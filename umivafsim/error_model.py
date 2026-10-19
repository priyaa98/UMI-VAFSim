import os
import sys
import json
import gzip
import hashlib
import time
import random
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pysam

from .utils import normalize_chrom, BASES

BASE_TO_INDEX = {b: i for i, b in enumerate(BASES)}

# --- ProcessPoolExecutor Worker State ---
_WORKER_BAMS = {}


def _init_worker(bam_paths, bam_threads):
    """Initializer hook called once when worker processes are created."""
    global _WORKER_BAMS
    _WORKER_BAMS = {p: pysam.AlignmentFile(p, "rb", threads=bam_threads) for p in bam_paths}


class BackgroundNoiseModel:
    """Empirical per-position substitution and quality profile learned from real control BAMs."""

    def __init__(self):
        self.position_probs = {}
        self.position_qual = {}
        self.cycle_qual_profile = {}
        self.cycle_qual_array = np.empty(0, dtype=np.float32)
        self.built_from = []
        self.provenance = None

    def _sync_cycle_array(self):
        """Converts dictionary profile into a dense 1D NumPy lookup array."""
        if not self.cycle_qual_profile:
            self.cycle_qual_array = np.empty(0, dtype=np.float32)
            return
        max_cycle = max(self.cycle_qual_profile.keys())
        arr = np.full(max_cycle + 1, -1.0, dtype=np.float32)
        for cyc, qval in self.cycle_qual_profile.items():
            arr[cyc] = qval
        self.cycle_qual_array = arr

    def to_dict(self):
        return {
            "position_probs": {f"{c}:{p}": v.tolist() for (c, p), v in self.position_probs.items()},
            "position_qual": {f"{c}:{p}": v for (c, p), v in self.position_qual.items()},
            "cycle_qual_profile": self.cycle_qual_profile,
            "built_from": self.built_from,
            "provenance": self.provenance,
        }

    @classmethod
    def from_dict(cls, d):
        m = cls()
        for k, v in d["position_probs"].items():
            c, p = k.rsplit(":", 1)
            m.position_probs[(c, int(p))] = np.array(v)
        for k, v in d["position_qual"].items():
            c, p = k.rsplit(":", 1)
            m.position_qual[(c, int(p))] = v
        m.cycle_qual_profile = {int(k): float(v) for k, v in d["cycle_qual_profile"].items()}
        # FIX (#2): Populate cycle_qual_array AFTER assigning cycle_qual_profile
        m._sync_cycle_array()
        m.built_from = d.get("built_from", [])
        m.provenance = d.get("provenance")
        return m

    def save(self, path):
        with gzip.open(path, "wt") as f:
            json.dump(self.to_dict(), f)

    @classmethod
    def load(cls, path):
        with gzip.open(path, "rt") as f:
            return cls.from_dict(json.load(f))


def compute_provenance(bam_paths, bed_regions=None):
    """Fingerprint of the control BAMs + BED regions a model is/was built from."""
    normalized_bams = sorted(os.path.abspath(p) for p in bam_paths)
    normalized_regions = sorted(list(r) for r in bed_regions) if bed_regions else []
    payload = json.dumps({"bams": normalized_bams, "regions": normalized_regions}, sort_keys=True)
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def merge_regions(regions):
    """Merges overlapping/adjacent (chrom, start, end) intervals per chrom."""
    by_chrom = defaultdict(list)
    for chrom, start, end in regions:
        by_chrom[chrom].append((start, end))

    merged = []
    for chrom, spans in by_chrom.items():
        if not spans:
            continue
        spans.sort()
        cur_start, cur_end = spans[0]
        for start, end in spans[1:]:
            if start <= cur_end:
                cur_end = max(cur_end, end)
            else:
                merged.append((chrom, cur_start, cur_end))
                cur_start, cur_end = start, end
        merged.append((chrom, cur_start, cur_end))

    return merged


def tile_regions(regions, window_size=100_000):
    """Splits large intervals into fixed-size windows to bound memory usage."""
    tiled = []
    for chrom, start, end in regions:
        for w_start in range(start, end, window_size):
            w_end = min(w_start + window_size, end)
            tiled.append((chrom, w_start, w_end))
    return tiled


def _median_from_histogram(hist):
    """Computes median Phred quality score directly from a count histogram array."""
    total = hist.sum()
    if total == 0:
        return 0.0
    cumulative = np.cumsum(hist)
    cutoff = total / 2.0
    return float(np.searchsorted(cumulative, cutoff))


def _derive_subseed(base_seed, label):
    """Derives a deterministic subseed integer from a base seed string/number and label."""
    payload = f"{base_seed}:{label}".encode("utf-8")
    return int(hashlib.sha256(payload).hexdigest()[:8], 16)


def _sample_cycle_worker(bam_path, units_chunk, reads_per_unit, max_reads_per_chunk):
    """Worker task to sample cycle-dependent qualities across a chunk of regions."""
    bam = _WORKER_BAMS[bam_path]
    local_hist = defaultdict(lambda: np.zeros(61, dtype=np.int32))
    sampled_count = 0

    for norm_chrom, start, end in units_chunk:
        if sampled_count >= max_reads_per_chunk:
            break

        try:
            fetch_iter = bam.fetch(norm_chrom, start, end)
        except ValueError:
            continue

        unit_sampled = 0
        for read in fetch_iter:
            if (read.is_unmapped or read.is_duplicate or read.is_qcfail or 
                    read.is_secondary or read.is_supplementary):
                continue

            quals = read.query_qualities
            if quals is None or len(quals) == 0:
                continue

            q_len = len(quals)
            is_rev = read.is_reverse
            clipped_quals = np.clip(quals, 0, 60)

            for qpos, bounded_qual in enumerate(clipped_quals):
                cycle = qpos if not is_rev else (q_len - 1 - qpos)
                local_hist[cycle][bounded_qual] += 1

            unit_sampled += 1
            sampled_count += 1

            if unit_sampled >= reads_per_unit:
                break

    return dict(local_hist)


def _sample_cycle_quals_parallel(executor, bam_paths, work_units, num_workers=4, max_reads_per_bam=50000, seed=42):
    """Parallelized stratified sampling that chunks regions across workers and batches intervals."""
    if not work_units:
        return {}

    with pysam.AlignmentFile(bam_paths[0], "rb") as first_bam:
        ref_lengths = dict(zip(first_bam.references, first_bam.lengths))

    coalesced_units = merge_regions([
        (
            chrom, 
            max(0, start - 100), 
            min(end + 100, ref_lengths.get(chrom, end + 100))
        ) 
        for chrom, start, end in work_units
    ])

    subseed = _derive_subseed(seed, "cycle_qual_sampling")
    rng = random.Random(subseed)
    units_shuffled = list(coalesced_units)
    rng.shuffle(units_shuffled)

    reads_per_unit = max(5, max_reads_per_bam // max(1, len(units_shuffled)))

    chunk_size = max(1, len(units_shuffled) // num_workers)
    region_chunks = [units_shuffled[i:i + chunk_size] for i in range(0, len(units_shuffled), chunk_size)]
    max_reads_per_chunk = max(1, max_reads_per_bam // len(region_chunks))

    futures = []
    for bam_path in bam_paths:
        for chunk in region_chunks:
            futures.append(
                executor.submit(_sample_cycle_worker, bam_path, chunk, reads_per_unit, max_reads_per_chunk)
            )

    combined_hist = defaultdict(lambda: np.zeros(61, dtype=np.int32))
    for future in as_completed(futures):
        local_hist = future.result()
        for cycle, hist in local_hist.items():
            combined_hist[cycle] += hist

    return {cyc: _median_from_histogram(hist) for cyc, hist in combined_hist.items()}


def _scan_chunk_worker(bam_path, chunk_units, min_mapq, min_baseq):
    """Worker function executing array-backed scanning over small windowed chunks."""
    FLAG_EXCLUDE = 0x4 | 0x100 | 0x200 | 0x400 | 0x800
    bam = _WORKER_BAMS[bam_path]

    sparse_results = {}

    for norm_chrom, start, end in chunk_units:
        length = end - start
        counts = np.zeros((length, 4), dtype=np.int64)
        quals_hist = np.zeros((length, 61), dtype=np.int32)

        try:
            pileup_columns = bam.pileup(
                norm_chrom, start, end,
                truncate=True,
                min_base_quality=min_baseq,
                min_mapping_quality=min_mapq,
                flag_filter=FLAG_EXCLUDE,
                ignore_overlaps=False,
            )
            
            for col in pileup_columns:
                offset = col.reference_pos - start
                if offset < 0 or offset >= length:
                    continue

                for pileup_read in col.pileups:
                    if pileup_read.is_del or pileup_read.is_refskip:
                        continue

                    qpos = pileup_read.query_position
                    if qpos is None:
                        continue

                    read = pileup_read.alignment
                    if read.is_unmapped or read.query_sequence is None or read.query_qualities is None:
                        continue

                    if qpos >= len(read.query_sequence) or qpos >= len(read.query_qualities):
                        continue

                    base = read.query_sequence[qpos].upper()
                    base_idx = BASE_TO_INDEX.get(base)
                    if base_idx is None:
                        continue

                    raw_q = read.query_qualities[qpos]
                    bounded_qual = 0 if raw_q < 0 else (60 if raw_q > 60 else raw_q)

                    counts[offset, base_idx] += 1
                    quals_hist[offset, bounded_qual] += 1

        except ValueError:
            pass

        covered_indices = np.flatnonzero(counts.sum(axis=1))
        for idx in covered_indices:
            key = (norm_chrom, start + idx)
            sparse_results[key] = (counts[idx], quals_hist[idx])

    return sparse_results


def build_background_model(bam_paths, bed_regions=None, min_mapq=20, min_baseq=10,
                            variant_vaf_threshold=0.01, min_supporting_reads=3,
                            num_workers=4, window_size=100_000, bam_threads=1, seed=42):
    bam_paths = [os.path.abspath(os.path.expanduser(p)) for p in bam_paths]
    """Parallelized, sparse-IPC BackgroundNoiseModel builder."""
    has_bed = bool(bed_regions)
    mode_str = f"across {len(bed_regions):,} region(s)" if has_bed else "across ALL reference contigs"
    print(f"[-] Building background noise model from {len(bam_paths)} control BAM(s) {mode_str}...")

    for path in bam_paths:
        if not os.path.exists(path):
            sys.exit(f"[FATAL] Control BAM file not found: {path}")

    if num_workers > 1 and bam_threads > 1:
        print(f"  [!] Note: Setting bam_threads=1 to prevent oversubscription with num_workers={num_workers}.")
        bam_threads = 1

    model = BackgroundNoiseModel()
    model.built_from = list(bam_paths)
    model.provenance = compute_provenance(bam_paths, bed_regions)

    first_bam = pysam.AlignmentFile(bam_paths[0], "rb", threads=bam_threads)
    bam_refs = set(first_bam.references)

    if has_bed:
        raw_units = []
        for chrom, start, end in bed_regions:
            norm_chrom = normalize_chrom(chrom, bam_refs)
            if norm_chrom in bam_refs:
                raw_units.append((norm_chrom, start, end))
        work_units = merge_regions(raw_units)
        if len(work_units) < len(raw_units):
            print(f"  [+] Merged {len(raw_units):,} BED intervals into {len(work_units):,} region(s).")
    else:
        work_units = [(ref, 0, length) for ref, length in zip(first_bam.references, first_bam.lengths)]
    first_bam.close()

    tiled_units = tile_regions(work_units, window_size=window_size)

    pos_counts = defaultdict(lambda: np.zeros(4, dtype=np.int64))
    pos_quals_hist = defaultdict(lambda: np.zeros(61, dtype=np.int32))

    chunk_size = max(1, len(tiled_units) // (num_workers * 4)) if tiled_units else 1
    chunks = [tiled_units[i:i + chunk_size] for i in range(0, len(tiled_units), chunk_size)]

    t_start = time.monotonic()
    total_chunks = len(chunks) * len(bam_paths)
    completed_chunks = 0

    with ProcessPoolExecutor(
        max_workers=num_workers, 
        initializer=_init_worker, 
        initargs=(bam_paths, bam_threads)
    ) as executor:

        print("  [+] Sampling cycle-dependent quality profile across regions (parallelized)...")
        model.cycle_qual_profile = _sample_cycle_quals_parallel(
            executor, bam_paths, work_units, num_workers=num_workers, seed=seed
        )
        model._sync_cycle_array()

        futures = []
        for bam_path in bam_paths:
            for chunk in chunks:
                futures.append(executor.submit(_scan_chunk_worker, bam_path, chunk, min_mapq, min_baseq))

        for future in as_completed(futures):
            sparse_chunk_res = future.result()
            for key, (c_arr, q_hist) in sparse_chunk_res.items():
                pos_counts[key] += c_arr
                pos_quals_hist[key] += q_hist

            completed_chunks += 1
            if completed_chunks % max(1, total_chunks // 5) == 0 or completed_chunks == total_chunks:
                elapsed = time.monotonic() - t_start
                print(f"  [+] Scanning progress: {completed_chunks}/{total_chunks} chunks completed ({elapsed:.1f}s elapsed)")

    print("  [+] Filtering variant candidates from noise profile...")
    for key, counts in pos_counts.items():
        depth = counts.sum()
        if depth == 0:
            continue

        major_idx = counts.argmax()
        non_major = depth - counts[major_idx]
        vaf = non_major / depth

        if vaf <= variant_vaf_threshold or non_major < min_supporting_reads:
            model.position_probs[key] = counts / depth
            model.position_qual[key] = _median_from_histogram(pos_quals_hist[key])

    print(f"[+] Background noise model built: {len(model.position_probs):,} genomic positions tracked, "
          f"{len(model.cycle_qual_profile):,} read-cycle quality bins.")

    return model
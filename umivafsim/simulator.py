import os
import sys
import gzip
import io
import hashlib
import shutil
import subprocess
import pysam
import random
import csv
import statistics
import numpy as np
from collections import defaultdict, OrderedDict
from contextlib import ExitStack
from .utils import (
    sample_negative_binomial_multiplier,
    calculate_injected_qual,
    calculate_injected_qual_v2,
    apply_sequencing_error_model,
    apply_sequencing_error_model_v2,
    pick_orientation,
    read_direction,
    strand_bias_reachable_range,
    normalize_chrom,
    parse_mi_duplex,
    parse_rx_duplex,
    split_rx_pair,
    extract_umi_sequence,
    build_fastq_header,
    reverse_complement,
    DUPLEX_DELIMITER,
    PHRED_TRANS_TABLE
)


def strip_boost_suffix(qname):
    """Removes a trailing _boost_N. Use for INTERNAL lookups (locus_reads, read_boost_map, mate pairing)."""
    return qname.rsplit("_boost_", 1)[0] if "_boost_" in qname else qname


def sanitize_qname(qname):
    """Canonical EMITTED read name: _boost_N stripped, truncated at the first '+'.
    Use for FASTQ names and the mutated-reads sidecar so the two always match 1:1.
    Do NOT use for internal dict lookups, which are keyed on strip_boost_suffix()."""
    clean = strip_boost_suffix(qname)
    plus_idx = clean.find("+")
    return clean[:plus_idx] if plus_idx != -1 else clean


# --- Mate pairing for the coordinate-sorted input BAM -------------------------------------------------------------
# Mates on the same contig are close together in a coordinate-sorted file, so a bounded in-memory buffer finds them.
# Mates on different contigs (or very far apart) are not: the first mate would wait in the buffer until it overflows
# and be flushed as an orphan before the second mate arrives. Such reads are spilled to a temporary BAM, name-sorted,
# and paired from there, so they stay proper pairs instead of becoming orphans.
MAX_PENDING_MATES = 100_000
DISTANT_MATE_BP = 1_000_000


def mate_is_distant(read):
    """True when the mate is too far away in a coordinate-sorted BAM to be held in memory until it arrives."""
    if not read.is_paired or read.mate_is_unmapped or read.next_reference_id < 0 or read.reference_id < 0:
        return False
    if read.next_reference_id != read.reference_id:
        return True
    return abs(read.next_reference_start - read.reference_start) > DISTANT_MATE_BP


def iter_mate_pairs(bam, spill_dir, threads=1, max_pending=MAX_PENDING_MATES, stats=None):
    """Yields ("pair", r1, r2) or ("orphan", read, None) for the primary reads of a coordinate-sorted BAM.

    Near mates are paired in memory (as before). Reads whose mate is distant are spilled to a temporary BAM
    in spill_dir, name-sorted, and paired afterwards. Only reads whose mate is really absent are orphans."""
    if stats is None:
        stats = {}
    stats.setdefault("distant_pairs_paired", 0)
    pending = OrderedDict()
    spill_path = None
    spill_out = None
    tmp_dir = None
    spilled = 0

    try:
        for read in bam.fetch(until_eof=True):
            if read.is_secondary or read.is_supplementary or not read.query_sequence:
                continue

            clean_qname = strip_boost_suffix(read.query_name)
            mate = pending.pop(clean_qname, None)

            if mate is not None:
                if mate.is_read1:
                    yield "pair", mate, read
                else:
                    yield "pair", read, mate
            elif mate_is_distant(read):
                if spill_out is None:
                    import tempfile
                    tmp_dir = tempfile.mkdtemp(prefix="spikein_mates_", dir=spill_dir)
                    spill_path = os.path.join(tmp_dir, "distant_mates.bam")
                    spill_out = pysam.AlignmentFile(spill_path, "wb", template=bam)
                spill_out.write(read)
                spilled += 1
            else:
                if len(pending) >= max_pending:
                    _, orphan_read = pending.popitem(last=False)
                    yield "orphan", orphan_read, None
                pending[clean_qname] = read

        while pending:
            _, orphan_read = pending.popitem(last=False)
            yield "orphan", orphan_read, None

        if spill_out is not None:
            spill_out.close()
            spill_out = None
            sorted_path = os.path.join(tmp_dir, "distant_mates.namesorted.bam")
            pysam.sort("-n", "-@", str(max(1, threads)), "-T", os.path.join(tmp_dir, "sort"), "-o", sorted_path, spill_path)
            os.remove(spill_path)

            with pysam.AlignmentFile(sorted_path, "rb") as sorted_bam:
                group_name, group = None, []
                for read in sorted_bam.fetch(until_eof=True):
                    name = strip_boost_suffix(read.query_name)
                    if group and name != group_name:
                        yield from _emit_name_group(group, stats)
                        group = []
                    group_name = name
                    group.append(read)
                if group:
                    yield from _emit_name_group(group, stats)
    finally:
        if spill_out is not None:
            spill_out.close()
        if tmp_dir is not None:
            shutil.rmtree(tmp_dir, ignore_errors=True)


def _emit_name_group(group, stats):
    """One name's reads from the name-sorted spill: one R1 + one R2 become a pair, anything else is an orphan."""
    r1s = [r for r in group if r.is_read1]
    r2s = [r for r in group if not r.is_read1]
    used = []
    if r1s and r2s:
        stats["distant_pairs_paired"] += 1
        used = [r1s[0], r2s[0]]
        yield "pair", r1s[0], r2s[0]
    for r in group:
        if not any(r is u for u in used):
            yield "orphan", r, None


def resolve_header_umi_token(read, umi_opts, duplex_opts):
    """Mate-independent UMI token for --umi-placement header; R1 and R2 get the identical string.
    duplex: the full 'u1-u2' RX pair. single: the UMI carried by the designated UMI read."""
    umi_len = umi_opts.get("umi_len", 12)
    if duplex_opts.get("mode") == "duplex":
        rx_val = read.get_tag("RX") if read.has_tag("RX") else None
        if rx_val and DUPLEX_DELIMITER in rx_val:
            u1, e1 = split_rx_pair(rx_val, 1, umi_len)
            u2, e2 = split_rx_pair(rx_val, 2, umi_len)
            if e1 is None and e2 is None:
                return f"{u1}{DUPLEX_DELIMITER}{u2}"
        return extract_umi_sequence(read, 1, umi_len)
    return extract_umi_sequence(read, umi_opts.get("umi_read") or 1, umi_len)


def single_mode_family_key(read, rx_val, mi_val):
    """Single-mode molecule key. MI is already a per-molecule ID, so prefer it. Otherwise RX alone would
    merge unrelated fragments that share a UMI, so scope it by contig + fragment (template) start.
    Both mates compute min(own start, mate start), so R1 and R2 land in the same family.
    Mates that are unmapped / on another contig fall back to their own start (single-read family)."""
    if mi_val:
        return str(mi_val)
    if read.next_reference_id == read.reference_id and read.next_reference_start >= 0:
        tmpl_start = min(read.reference_start, read.next_reference_start)
    else:
        tmpl_start = read.reference_start
    return f"{rx_val}@{read.reference_id}:{tmpl_start}"


class FastqChunkWriter:
    """Buffers byte outputs into chunks for parallel pigz subprocess or standard file streams.

    Use as a context manager: a clean exit calls close() (flush + verify pigz exit code, raises
    RuntimeError on failure); an exit via exception calls abort() (kill pigz, delete partial file,
    never raises) so the original exception is never masked."""

    def __init__(self, target_path, threads=4, batch_size=50000):
        self.target_path = target_path
        self.batch_size = batch_size
        self.buffer = []
        self._closed = False
        self.proc = None
        self.out_fh = None
        self.stream = None

        is_gz = target_path.endswith(".gz")
        pigz_bin = shutil.which("pigz") if is_gz else None

        if is_gz and pigz_bin:
            cmd = [pigz_bin, "-p", str(threads), "-1"]
            self.out_fh = open(target_path, "wb")
            self.proc = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=self.out_fh, bufsize=8 * 1024 * 1024)
            self.stream = self.proc.stdin
            self.mode = "subprocess_pigz"
        elif is_gz:
            print(f"[WARNING] 'pigz' executable not found on PATH. Falling back to single-threaded gzip for '{target_path}'. "
                  f"Performance will be degraded. Install pigz for parallel compression.")
            self.stream = gzip.open(target_path, "wb", compresslevel=1)
            self.mode = "gzip_wrapper"
        else:
            self.stream = open(target_path, "wb", buffering=8 * 1024 * 1024)
            self.mode = "text_file"

    def write_record(self, record_bytes: bytes):
        self.buffer.append(record_bytes)
        if len(self.buffer) >= self.batch_size:
            self.flush()

    def flush(self):
        if not self.buffer:
            return
        payload = b"".join(self.buffer)
        self.buffer.clear()
        self.stream.write(payload)

    def _terminate(self):
        """Best-effort teardown of the subprocess and file handles. Never raises."""
        try:
            if self.proc is not None and self.proc.poll() is None:
                self.proc.kill()
        except Exception:
            pass
        for closer in (
            lambda: self.stream.close() if self.stream is not None else None,
            lambda: self.proc.wait() if self.proc is not None else None,
            lambda: self.out_fh.close() if self.out_fh is not None else None,
        ):
            try:
                closer()
            except Exception:
                pass

    def close(self):
        """Finalize output. Raises RuntimeError (not SystemExit) if the stream or pigz fails."""
        if self._closed:
            return
        self._closed = True
        try:
            self.flush()
            self.stream.close()
        except Exception as e:
            self._terminate()
            raise RuntimeError(f"Failed while finalizing '{self.target_path}': {e}") from e

        if self.mode == "subprocess_pigz":
            rc = self.proc.wait()
            self.out_fh.close()
            if rc != 0:
                raise RuntimeError(
                    f"pigz exited with code {rc} while writing '{self.target_path}'. "
                    f"Output file may be incomplete or corrupted."
                )

    def abort(self):
        """Failure path: drop buffered data, stop pigz, delete the partial file. Never raises."""
        self._closed = True
        try:
            self.buffer.clear()
        except Exception:
            pass
        self._terminate()
        try:
            if os.path.exists(self.target_path):
                os.remove(self.target_path)
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            self.close()
        else:
            self.abort()
        return False


class UmiEmbedTracker:
    def __init__(self):
        self.total_reads = 0
        self.failed_embeds = 0

    def log_result(self, is_fallback: bool):
        self.total_reads += 1
        if is_fallback:
            self.failed_embeds += 1

    def summarize(self):
        if self.total_reads > 0 and self.failed_embeds > 0:
            pct = (self.failed_embeds / self.total_reads) * 100
            print(f"[WARNING] UMI Embedding Fallbacks: {self.failed_embeds:,} / {self.total_reads:,} "
                  f"({pct:.2f}%) reads fell back to placeholder N-bases.")


def validate_bam_sort_and_index(bam_path, sample_size=10000):
    """Pre-flight check: confirms index existence, coordinate sort status, and paired flag validity."""
    print(f"[-] Pre-flight check: Validating sort order and index integrity for '{bam_path}'...")
    if not os.path.exists(bam_path):
        sys.exit(f"[FATAL] Alignment BAM file does not exist: {bam_path}")

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")
    except Exception as e:
        sys.exit(f"[FATAL] Failed to open BAM file '{bam_path}': {e}")

    if not bam.has_index():
        bam.close()
        sys.exit(f"[FATAL] Missing alignment index (.bai or .csi) for '{bam_path}'. Index required for region fetching.")

    header = bam.header.to_dict()
    hd = header.get("HD", {})
    so = hd.get("SO", "unknown")
    if so != "coordinate":
        print(f"[WARNING] BAM header SO tag is '{so}', expected 'coordinate'. Verifying order empirically...")

    sampled = 0
    paired_cnt = 0
    prev_chrom = None
    prev_pos = -1
    order_violations = 0

    for read in bam.fetch(until_eof=True):
        sampled += 1
        if read.is_paired:
            paired_cnt += 1

        if not read.is_unmapped:
            curr_chrom = read.reference_id
            curr_pos = read.reference_start
            if prev_chrom is not None:
                if curr_chrom < prev_chrom or (curr_chrom == prev_chrom and curr_pos < prev_pos):
                    order_violations += 1
            prev_chrom = curr_chrom
            prev_pos = curr_pos

        if sampled >= sample_size:
            break

    bam.close()

    if sampled == 0:
        sys.exit(f"[FATAL] Alignment file '{bam_path}' contains zero reads.")

    if (paired_cnt / sampled) < 0.5:
        sys.exit(f"[FATAL] BAM file reads lack paired-end flag ({paired_cnt}/{sampled} paired). Dual FASTQ pipeline expects paired data.")

    if order_violations > 0:
        sys.exit(f"[FATAL] Input BAM file is not coordinate-sorted ({order_violations} out-of-order alignments detected in sample).")

    print("[+] Pre-flight BAM sort order, index, and paired-end checks passed.")


def assert_umi_tagged(bam_path, resolved_variants, probe_window=500, mode="single", umi_len=12, max_sampled_reads=500):
    """Validates presence and integrity of UMI tags using bounded sampling across variant loci."""
    bam = pysam.AlignmentFile(bam_path, "rb")
    bam_refs = set(bam.references)

    total_sampled = 0
    rx_present_count = 0
    rx_hyphen_count = 0
    rx_valid_len_count = 0
    mi_suffixed_count = 0
    embeddable_umi_count = 0

    # No variants (--no-spike): sample the first reads of the BAM instead of reads around variant loci.
    for v in (resolved_variants if resolved_variants else [None]):
        if total_sampled >= max_sampled_reads:
            break

        if v is None:
            try:
                fetch_iter = bam.fetch()
            except ValueError:
                fetch_iter = []
        else:
            norm_chrom = normalize_chrom(v["chrom"], bam_refs)
            if norm_chrom not in bam_refs:
                bam.close()
                sys.exit(f"[FATAL] Target chromosome '{v['chrom']}' (normalized: '{norm_chrom}') is not present in BAM header.")

            pos = v["pos"]
            try:
                fetch_iter = bam.fetch(norm_chrom, max(0, pos - probe_window), pos + probe_window)
            except ValueError:
                fetch_iter = []

        for read in fetch_iter:
            if read.is_secondary or read.is_supplementary:
                continue
            
            total_sampled += 1

            rx_val = read.get_tag("RX") if read.has_tag("RX") else None
            mi_val = read.get_tag("MI") if read.has_tag("MI") else None

            extracted_umi = extract_umi_sequence(read, 1, umi_len)
            if extracted_umi:
                embeddable_umi_count += 1

            if rx_val:
                rx_present_count += 1
                if DUPLEX_DELIMITER in rx_val:
                    rx_hyphen_count += 1
                    _, err = split_rx_pair(rx_val, 1, umi_len)
                    if err is None:
                        rx_valid_len_count += 1

            if mi_val:
                _, strand = parse_mi_duplex(mi_val)
                if strand is not None:
                    mi_suffixed_count += 1

            if total_sampled >= max_sampled_reads:
                break

    bam.close()

    if total_sampled == 0:
        sys.exit("[FATAL] No reads found across target loci (or in the BAM, with --no-spike) to validate UMI tagging.")

    if mode == "single":
        if (embeddable_umi_count / total_sampled) < 0.90:
            sys.exit(f"[FATAL] Missing UMI tags. Only {embeddable_umi_count}/{total_sampled} sampled reads carry extractable UMI tags.")
        print(f"[+] Upfront UMI Precondition Validated (SINGLE mode): {embeddable_umi_count}/{total_sampled} tagged & embeddable.")
        return "SINGLE"

    has_hyphen_structure = (rx_hyphen_count / total_sampled) >= 0.90
    has_mi_structure = (mi_suffixed_count / total_sampled) >= 0.90

    if not has_hyphen_structure and not has_mi_structure:
        sys.exit("[FATAL] Input data is not duplex-tagged (Neither hyphenated RX nor suffixed MI detected across >=90% of sampled reads).")

    if has_hyphen_structure and (rx_valid_len_count / total_sampled) < 0.90:
        sys.exit(f"[FATAL] RX UMI barcode length mismatch. Valid hyphenated RX tag ratio ({rx_valid_len_count}/{total_sampled}) is below 90%.")

    if (embeddable_umi_count / total_sampled) < 0.90:
        sys.exit(f"[FATAL] Tool requires embeddable UMI sequences. Only {embeddable_umi_count}/{total_sampled} sampled reads carry extractable UMIs.")

    strand_source = "MI" if has_mi_structure else "RX_pseudo"
    if strand_source == "RX_pseudo":
        print("\n[WARNING] ATTENTION: No physical strand suffix found in MI tags. Falling back to RX lexicographical tie-breaking.\n")

    print(f"[+] Upfront UMI Precondition Validated (DUPLEX mode, strand_source='{strand_source}'): {embeddable_umi_count}/{total_sampled} tagged & embeddable.")
    return strand_source


def check_bam_integrity_and_variants(bam_path, chrom, pos, ref_base, force_noise=False, alt_base=None, noise_threshold=0.01, bg=None):
    if not os.path.exists(bam_path):
        sys.exit(f"[FATAL] Coordinate BAM file does not exist at path: {bam_path}")

    bam = pysam.AlignmentFile(bam_path, "rb")
    bam_refs = set(bam.references)
    norm_chrom = normalize_chrom(chrom, bam_refs)

    if norm_chrom not in bam_refs:
        bam.close()
        sys.exit(f"[FATAL] Target chromosome '{chrom}' (normalized: '{norm_chrom}') is not present in BAM header.")

    pre_existing_variants, total_overlapping_bases, pre_existing_alt = 0, 0, 0
    alt_upper = alt_base.upper() if alt_base else None
    target_ref_pos = pos - 1
    
    try:
        fetch_iter = bam.fetch(norm_chrom, target_ref_pos, pos)
    except ValueError:
        fetch_iter = []

    for read in fetch_iter:
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        seq = read.query_sequence
        if not seq:
            continue
        for read_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if ref_pos == target_ref_pos and read_pos is not None and 0 <= read_pos < len(seq):
                total_overlapping_bases += 1
                if seq[read_pos].upper() != ref_base.upper():
                    pre_existing_variants += 1
                if alt_upper is not None and seq[read_pos].upper() == alt_upper:
                    pre_existing_alt += 1
    bam.close()

    # Background of the UNSPIKED reads at this site, read level (reported in the ledger; does not change the gate below).
    if bg is not None:
        bg.update({"total_reads": total_overlapping_bases, "nonref_reads": pre_existing_variants,
                   "alt_reads": pre_existing_alt})

    if total_overlapping_bases == 0:
        print(f"[WARNING] Zero sequence coverage found at {norm_chrom}:{pos}.")
        return 0

    variant_ratio = pre_existing_variants / total_overlapping_bases
    if variant_ratio > noise_threshold:
        print(f"[WARNING] Pre-existing noise ({variant_ratio*100:.2f}%) detected at {norm_chrom}:{pos}!")
        if not force_noise:
            print(f"[WARNING] Excluding locus {norm_chrom}:{pos} from processing due to background noise > {noise_threshold*100:.1f}%.")
            return -1

    return total_overlapping_bases


def scan_target_families_multi(bam_path, variants_list, boost_opts, duplex_opts, stochastic=True, global_seed=42):
    """Pass 1: Scans molecular families and assigns mutant status across multiple variants concurrently."""
    print(f"[-] Pass 1: Scanning target molecular families for {len(variants_list)} variant(s)...")

    locus_reads = defaultdict(dict)
    read_boost_map = {}
    variant_stats = {}

    bam = pysam.AlignmentFile(bam_path, "rb")
    bam_refs = set(bam.references)

    for idx, v_info in enumerate(variants_list):
        vid = v_info["cosmic_id"]
        if v_info.get("excluded_due_to_noise", False):
            variant_stats[vid] = {
                "population_size": 0, "chosen_families": 0, "realized_vaf": 0.0,
                "complete_duplex": 0, "single_strand_duplex": 0, "unparseable_tags": 0,
                "status": "EXCLUDED_NOISE", "selection_rule": "N/A"
            }
            continue

        chrom = normalize_chrom(v_info["chrom"], bam_refs)
        pos = v_info["pos"]
        vaf = v_info["vaf"]
        spike_mode = v_info["spike_mode"]

        seed_bytes = f"{global_seed}:{vid}:{idx}".encode("utf-8")
        var_subseed = int.from_bytes(hashlib.sha256(seed_bytes).digest()[:4], byteorder="big")
        rng = np.random.RandomState(var_subseed)
        py_rng = random.Random(var_subseed)

        molecule_groups = defaultdict(lambda: defaultdict(set)) if duplex_opts["mode"] == "duplex" else defaultdict(set)
        read_strand_map = {}
        # fam_id -> [(clean_qname, strand_key, is_reverse)], one entry per READ RECORD covering the site
        # (a pair whose mates both cover it contributes two entries). Feeds per-orientation planning.
        family_reads = defaultdict(list)
        unparseable_tags_count = 0
        target_ref_pos = pos - 1
        alt_upper = str(v_info.get("alt", "")).upper()
        # Background per family: [reads covering the site, reads already carrying the alt]. A family counts as carrying the
        # background alt only when at least half of its reads do (a real molecule-level variant); a single erroneous read
        # in a family of ~50 does not, because consensus calling would outvote it.
        bg_fam_counts = defaultdict(lambda: [0, 0])
        # Native read pair type per family at this site: fam_id -> [reads that are F1R2, reads that are F2R1].
        # F1R2 = read 1 forward / read 2 reverse; F2R1 = read 1 reverse / read 2 forward. Reads that are not part of a
        # pair carry no pair type and are not counted. Descriptive only: never used for family selection.
        fam_pt_reads = defaultdict(lambda: [0, 0])
        # Distance (in reference bases, within the aligned block) from the site to the read's 5' end and 3' end, per
        # covering read: fam_id -> [(strand_key, read_num, d5, d3)], read_num = 1, 2 or 0 (not part of a pair). The 5'
        # end is the aligned start of a forward read and the aligned end of a reverse read. Descriptive only; shows how
        # close spiked bases sit to the read ends (primers). Only d5 is reported so far; d3 is kept for later use.
        fam_read_dist = defaultdict(list)

        try:
            fetch_iter = bam.fetch(chrom, target_ref_pos, pos)
        except ValueError:
            fetch_iter = []

        for read in fetch_iter:
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue

            qname = read.query_name
            clean_qname = strip_boost_suffix(qname)

            rx_val = read.get_tag("RX") if read.has_tag("RX") else None
            mi_val = read.get_tag("MI") if read.has_tag("MI") else None

            if not rx_val and not mi_val:
                continue

            seq = read.query_sequence
            if not seq:
                continue

            hit_pos = next((rp for rp, rf in read.get_aligned_pairs(matches_only=True)
                            if rf == target_ref_pos and rp is not None and 0 <= rp < len(seq)), None)
            if hit_pos is not None:
                read_has_alt = bool(alt_upper) and seq[hit_pos].upper() == alt_upper
                if duplex_opts["mode"] == "duplex":
                    if duplex_opts["strand_source"] == "MI":
                        fam_id, strand_key = parse_mi_duplex(mi_val)
                    else:
                        fam_id, strand_key = parse_rx_duplex(rx_val)

                    if fam_id is None or strand_key is None:
                        unparseable_tags_count += 1
                    else:
                        molecule_groups[fam_id][strand_key].add(clean_qname)
                        read_strand_map[clean_qname] = (fam_id, strand_key)
                        family_reads[fam_id].append((clean_qname, strand_key, read.is_reverse))
                        bg_fam_counts[fam_id][0] += 1
                        bg_fam_counts[fam_id][1] += int(read_has_alt)
                        if read.is_read1 or read.is_read2:
                            fam_pt_reads[fam_id][0 if (read.is_read1 != read.is_reverse) else 1] += 1
                        _ds, _de = target_ref_pos - read.reference_start, read.reference_end - 1 - target_ref_pos
                        _rn = 1 if read.is_read1 else (2 if read.is_read2 else 0)
                        fam_read_dist[fam_id].append((strand_key, _rn, _de, _ds) if read.is_reverse else (strand_key, _rn, _ds, _de))
                else:
                    umi_key = single_mode_family_key(read, rx_val, mi_val)
                    molecule_groups[umi_key].add(clean_qname)
                    read_strand_map[clean_qname] = (umi_key, "A")
                    family_reads[umi_key].append((clean_qname, "A", read.is_reverse))
                    bg_fam_counts[umi_key][0] += 1
                    bg_fam_counts[umi_key][1] += int(read_has_alt)
                    if read.is_read1 or read.is_read2:
                        fam_pt_reads[umi_key][0 if (read.is_read1 != read.is_reverse) else 1] += 1
                    _ds, _de = target_ref_pos - read.reference_start, read.reference_end - 1 - target_ref_pos
                    _rn = 1 if read.is_read1 else (2 if read.is_read2 else 0)
                    fam_read_dist[umi_key].append(("A", _rn, _de, _ds) if read.is_reverse else ("A", _rn, _ds, _de))

        raw_total_molecules = len(molecule_groups)
        complete_duplex_cnt = 0
        single_strand_duplex_cnt = 0
        eligible_fam_ids = []

        if duplex_opts["mode"] == "duplex":
            for fam_id, strands in molecule_groups.items():
                has_a, has_b = bool(strands.get("A")), bool(strands.get("B"))
                if has_a and has_b:
                    complete_duplex_cnt += 1
                    eligible_fam_ids.append(fam_id)
                elif has_a or has_b:
                    single_strand_duplex_cnt += 1
                    if spike_mode == "single-strand":
                        eligible_fam_ids.append(fam_id)

            sampling_pool = eligible_fam_ids
            population_size = len(eligible_fam_ids)
        else:
            sampling_pool = list(molecule_groups.keys())
            population_size = raw_total_molecules

        if population_size == 0 or len(sampling_pool) == 0:
            print(f"[WARNING] Variant {vid} at {chrom}:{pos} yielded 0 eligible molecular families. Skipping mutation.")
            variant_stats[vid] = {
                "population_size": 0, "chosen_families": 0, "realized_vaf": 0.0,
                "complete_duplex": 0, "single_strand_duplex": 0, "unparseable_tags": unparseable_tags_count,
                "status": "EXCLUDED_ZERO_COVERAGE", "selection_rule": "none"
            }
            continue

        selection_rule = "binomial"
        if stochastic:
            num_families = rng.binomial(n=population_size, p=vaf)
            if num_families == 0 and vaf > 0:
                num_families = 1
                selection_rule = "forced_min_one"
        else:
            num_families = 0 if vaf == 0.0 else max(1, round(population_size * vaf))
            selection_rule = "deterministic"

        num_families = min(num_families, len(sampling_pool))
        chosen_fam_ids = py_rng.sample(sampling_pool, num_families)

        # Strand to mutate per chosen family (duplex single-strand only). NOTE: py_rng consumption here is
        # unchanged from earlier versions, so the base family set is identical with and without strand bias.
        strand_pick = {}
        for fam_id in chosen_fam_ids:
            if duplex_opts["mode"] == "duplex" and spike_mode == "single-strand":
                strands = molecule_groups[fam_id]
                available_strands = [st for st in ["A", "B"] if strands.get(st)]
                strand_pick[fam_id] = py_rng.choice(available_strands) if available_strands else "A"
            else:
                strand_pick[fam_id] = None

        # Orientation bias = which ORIGINAL strand carries the alt (duplex single-strand spikes only): the requested
        # fraction of alt-carrying families whose alt sits on strand A (the rest on strand B). Families covered by one
        # strand only have no choice; the others are assigned with a separate rng, so the base family set is unchanged.
        ob_enabled = bool(v_info.get("orientation_bias_enabled", False))
        ob_info = {}
        if ob_enabled and duplex_opts["mode"] != "duplex":
            sys.exit(f"[FATAL] Variant {vid}: orientation_bias is duplex-only (each single-UMI family is one strand).")
        if ob_enabled:
            ob_target = float(v_info["orientation_bias"])
            ob_rng = random.Random(int.from_bytes(hashlib.sha256(f"{global_seed}:{vid}:ob_strand".encode()).digest()[:4], "big"))
            forced_a, forced_b, free_ob = [], [], []
            for fam_id in chosen_fam_ids:
                strands = molecule_groups[fam_id]
                has_a, has_b = bool(strands.get("A")), bool(strands.get("B"))
                if has_a and has_b:
                    free_ob.append(fam_id)
                elif has_b:
                    forced_b.append(fam_id)
                else:
                    forced_a.append(fam_id)
            n_ch = len(chosen_fam_ids)
            ob_lo = len(forced_a) / n_ch if n_ch else None
            ob_hi = (len(forced_a) + len(free_ob)) / n_ch if n_ch else None
            want_a = int(round(ob_target * n_ch))
            n_free_a = max(0, min(len(free_ob), want_a - len(forced_a)))
            shuffled_ob = list(free_ob)
            ob_rng.shuffle(shuffled_ob)
            for fam_id in forced_a:
                strand_pick[fam_id] = "A"
            for fam_id in forced_b:
                strand_pick[fam_id] = "B"
            for i, fam_id in enumerate(shuffled_ob):
                strand_pick[fam_id] = "A" if i < n_free_a else "B"
            planned_a = len(forced_a) + n_free_a
            ob_planned = (planned_a / n_ch) if n_ch else None
            if n_ch == 0:
                ob_status = "NO_ALT_FAMILIES"
            elif ob_target < ob_lo - 1e-9 or ob_target > ob_hi + 1e-9:
                ob_status = "UNREACHABLE"
            elif abs(ob_planned - ob_target) <= max(0.02, 1.0 / n_ch):
                ob_status = "MET"
            else:
                ob_status = "APPROX"
            ob_info = {
                "ob_enabled": True, "ob_requested": ob_target, "ob_min": ob_lo, "ob_max": ob_hi,
                "ob_planned_a": ob_planned, "ob_status": ob_status,
                "ob_families_a": planned_a, "ob_families_b": n_ch - planned_a,
                "ob_forced_a": len(forced_a), "ob_forced_b": len(forced_b), "ob_free": len(free_ob),
            }

        def cand_counts(fam_id, pick):
            """(forward, reverse) consensus-unit counts for this family: 1 if a consensus read of that aligned
            orientation can be built at the site, else 0. Steering works on consensus units (one per family and
            orientation) because the caller sees consensus reads, not raw reads. In duplex mode a unit needs reads
            from both strands, otherwise duplex consensus drops it."""
            by_orient = {"F": set(), "R": set()}
            for _qn, st_key, is_rev in family_reads.get(fam_id, ()):
                if pick is not None and st_key != pick:
                    continue
                by_orient["R" if is_rev else "F"].add(st_key)
            need = 2 if (duplex_opts["mode"] == "duplex" and pick is None) else 1
            return (1 if len(by_orient["F"]) >= need else 0), (1 if len(by_orient["R"]) >= need else 0)

        # Orientation composition of the site, recorded in every run (also when no strand bias is requested). Descriptive
        # only: it uses no random numbers and does not influence which families are chosen or mutated.
        # Pool = all eligible families (all strands); chosen = the families that carry the alt.
        pool_f_only = pool_r_only = pool_both = 0
        for f in sampling_pool:
            c_f, c_r = cand_counts(f, None)
            if c_r == 0 and c_f > 0:
                pool_f_only += 1
            elif c_f == 0 and c_r > 0:
                pool_r_only += 1
            elif c_f > 0 and c_r > 0:
                pool_both += 1
        chosen_f_only = chosen_r_only = chosen_both = 0
        for f in chosen_fam_ids:
            c_f, c_r = cand_counts(f, strand_pick[f])
            if c_r == 0 and c_f > 0:
                chosen_f_only += 1
            elif c_f == 0 and c_r > 0:
                chosen_r_only += 1
            elif c_f > 0 and c_r > 0:
                chosen_both += 1

        def _pt_class(fams):
            """Families whose reads at the site are all F1R2, all F2R1, or of both pair types."""
            only1 = only2 = both_pt = 0
            n1 = n2 = 0
            for f in fams:
                a, b = fam_pt_reads.get(f, (0, 0))
                n1 += a
                n2 += b
                if a and b:
                    both_pt += 1
                elif a:
                    only1 += 1
                elif b:
                    only2 += 1
            return only1, only2, both_pt, n1, n2

        pool_pt = _pt_class(sampling_pool)
        chosen_pt = _pt_class(chosen_fam_ids)
        all_n1 = sum(v[0] for v in fam_pt_reads.values())
        all_n2 = sum(v[1] for v in fam_pt_reads.values())
        def _dist_stats(fams, picks=None):
            """Distance from the site to the read's 5' end (aligned bases), per read number (1, 2): min and median
            over the covering reads of the given families (restricted to the picked strand if given). Descriptive only;
            no threshold is applied here."""
            d = {1: [], 2: [], "all": []}
            for f in fams:
                pk = picks.get(f) if picks else None
                for st_key, rn, d5, d3 in fam_read_dist.get(f, ()):
                    if rn in (1, 2) and (pk is None or st_key == pk):
                        d[rn].append(d5)
                        d["all"].append(d5)
            return {rn: ((min(v), statistics.median(v)) if v else None) for rn, v in d.items()}

        chosen_ds = _dist_stats(chosen_fam_ids, strand_pick)
        native_ds = _dist_stats(list(fam_read_dist.keys()))
        comp_info = {
            "sb_pool_fwd_only": pool_f_only, "sb_pool_rev_only": pool_r_only, "sb_pool_both": pool_both,
            "chosen_fwd_only": chosen_f_only, "chosen_rev_only": chosen_r_only, "chosen_both": chosen_both,
            "pool_f1r2_fams": pool_pt[0], "pool_f2r1_fams": pool_pt[1], "pool_both_pt_fams": pool_pt[2],
            "chosen_f1r2_fams": chosen_pt[0], "chosen_f2r1_fams": chosen_pt[1], "chosen_both_pt_fams": chosen_pt[2],
            "native_f1r2_reads": all_n1, "native_f2r1_reads": all_n2,
            "chosen_f1r2_reads": chosen_pt[3], "chosen_f2r1_reads": chosen_pt[4],
            "chosen_dist": chosen_ds, "native_dist": native_ds,
        }

        # Read-orientation strand bias (the alt is placed on forward- or reverse-aligned reads of a family).
        sb_enabled = bool(v_info.get("strand_bias_enabled", False))
        sb_target = float(v_info.get("fwd_bias_ratio", 0.5))
        fam_orient = {fam_id: None for fam_id in chosen_fam_ids}   # None = mutate every covering read
        sb_info = {}

        if sb_enabled:
            orient_rng = random.Random(int.from_bytes(hashlib.sha256(f"{global_seed}:{vid}:sb_orient".encode()).digest()[:4], "big"))

            # The chosen families are NOT changed: the VAF is defined by how many families carry the alt. Strand bias
            # only decides, inside those families, which consensus units (family x aligned orientation) keep the alt.
            # Reads of a unit that does not carry the alt stay in the data as reference; nothing is removed.
            counts_by_fam = {f: cand_counts(f, strand_pick[f]) for f in chosen_fam_ids}
            forced_f = [f for f in chosen_fam_ids if counts_by_fam[f] == (1, 0)]
            forced_r = [f for f in chosen_fam_ids if counts_by_fam[f] == (0, 1)]
            free = [f for f in chosen_fam_ids if counts_by_fam[f] == (1, 1)]
            n_f, n_r, n_b = len(forced_f), len(forced_r), len(free)
            for f in forced_f:
                fam_orient[f] = "F"
            for f in forced_r:
                fam_orient[f] = "R"

            # Free families start with the alt on both orientations (no skew). To move the forward fraction of alt
            # units to the requested value, the alt is withheld from one orientation of just enough free families.
            fu, ru = n_f + n_b, n_r + n_b
            n_drop, drop_mode = 0, None
            if fu + ru > 0:
                if sb_target > fu / (fu + ru):           # too few forward: withhold the alt from reverse units
                    n_drop, drop_mode = ru - fu * (1.0 - sb_target) / sb_target, "F"
                elif sb_target < fu / (fu + ru):         # too few reverse: withhold the alt from forward units
                    n_drop, drop_mode = fu - ru * sb_target / (1.0 - sb_target), "R"
                n_drop = max(0, min(n_b, int(round(n_drop))))
            shuffled = list(free)
            orient_rng.shuffle(shuffled)
            for f in shuffled[:n_drop]:
                fam_orient[f] = drop_mode            # keeps the alt on this orientation only
            alt_f = n_f + n_b - (n_drop if drop_mode == "R" else 0)
            alt_r = n_r + n_b - (n_drop if drop_mode == "F" else 0)
            alt_total = alt_f + alt_r
            potential_units = n_f + n_r + 2 * n_b

            lo, hi = strand_bias_reachable_range(counts_by_fam.values())
            planned = (alt_f / alt_total) if alt_total > 0 else None
            tol_eps = 1e-9
            if lo is None or hi is None or alt_total == 0:
                plan_status = "NO_ALT_READS"
            elif sb_target < lo - tol_eps or sb_target > hi + tol_eps:
                plan_status = "UNREACHABLE"
            elif abs(planned - sb_target) <= max(0.02, 1.0 / alt_total):
                plan_status = "MET"
            else:
                plan_status = "APPROX"

            sb_info = {
                "sb_enabled": True, "sb_requested": sb_target,
                "sb_min": lo, "sb_max": hi, "sb_planned_fwd": planned,
                "sb_plan_status": plan_status,
                "sb_alt_units_potential": potential_units, "sb_alt_units_planned": alt_total,
            }

        chosen_fam_ids_set = set(chosen_fam_ids)

        # qname -> orientation that carries the alt (None = every covering read). Boosted copies inherit it.
        mutant_orient = {}
        for fam_id in chosen_fam_ids:
            for qn, st_key, _is_rev in family_reads.get(fam_id, ()):
                pick = strand_pick.get(fam_id)
                if pick is not None and st_key != pick:
                    continue
                mutant_orient[qn] = fam_orient[fam_id]

        # Expected consensus-level VAF. After duplex/strand consensus each family yields one consensus read per
        # aligned orientation (F/R) covering the site, merging the reads of both strands. A unit carries the alt only
        # if every strand present in it carries it; a unit with mixed strands is discordant (masked by duplex consensus).
        cons_units = {}
        for _fam, _recs in family_reads.items():
            for _qn, _st_key, _is_rev in _recs:
                _o = read_direction(_is_rev)
                _is_alt = _qn in mutant_orient and mutant_orient[_qn] in (None, _o)
                _u = cons_units.setdefault((_fam, _o), {})
                _u[_st_key] = _u.get(_st_key, False) or _is_alt
        if duplex_opts["mode"] == "duplex":
            # duplex consensus (min-reads 1 1 1) only emits a read where both strands contribute
            cons_units = {k: u for k, u in cons_units.items() if len(u) >= 2}
        cons_total = len(cons_units)
        cons_alt = sum(1 for u in cons_units.values() if all(u.values()))
        cons_disc = sum(1 for u in cons_units.values() if any(u.values()) and not all(u.values()))
        cons_alt_f = sum(1 for (_f, _o), u in cons_units.items() if _o == "F" and all(u.values()))
        cons_info = {
            "cons_units_total": cons_total, "cons_units_alt": cons_alt, "cons_units_discordant": cons_disc,
            "cons_alt_fwd_fraction": (cons_alt_f / cons_alt) if cons_alt else None,
            "expected_consensus_vaf": (cons_alt / cons_total * 100.0) if cons_total else 0.0,
        }

        if boost_opts["enabled"]:
            if duplex_opts["mode"] == "duplex":
                for fam_id, strands in molecule_groups.items():
                    is_mutant_fam = fam_id in chosen_fam_ids_set
                    should_boost = (boost_opts["mode"] == "all-locus" or is_mutant_fam)

                    if duplex_opts["boost_linkage"] == "joint":
                        mult = sample_negative_binomial_multiplier(boost_opts["factor"], boost_opts["dispersion"], rng=rng) if should_boost else 1
                        for st_qnames in strands.values():
                            for qn in st_qnames:
                                read_boost_map[qn] = max(read_boost_map.get(qn, 1), mult)
                    else:
                        for st_qnames in strands.values():
                            mult = sample_negative_binomial_multiplier(boost_opts["factor"], boost_opts["dispersion"], rng=rng) if should_boost else 1
                            for qn in st_qnames:
                                read_boost_map[qn] = max(read_boost_map.get(qn, 1), mult)
            else:
                for tag, qnames in molecule_groups.items():
                    mult = sample_negative_binomial_multiplier(boost_opts["factor"], boost_opts["dispersion"], rng=rng) \
                        if (boost_opts["mode"] == "all-locus" or tag in chosen_fam_ids_set) else 1
                    for qn in qnames:
                        read_boost_map[qn] = max(read_boost_map.get(qn, 1), mult)

        all_locus_qnames = set()
        if duplex_opts["mode"] == "duplex":
            for strands in molecule_groups.values():
                for st_qnames in strands.values():
                    all_locus_qnames.update(st_qnames)
        else:
            for qnames in molecule_groups.values():
                all_locus_qnames.update(qnames)

        for qn in all_locus_qnames:
            fam_id, strand_key = read_strand_map.get(qn, (qn, "A"))
            locus_reads[qn][vid] = {
                "is_mutant": qn in mutant_orient,
                "mut_orient": mutant_orient.get(qn),
                "fam_id": fam_id,
                "strand_key": strand_key
            }

        realized_vaf = (len(chosen_fam_ids) / population_size * 100.0) if population_size > 0 else 0.0
        variant_stats[vid] = {
            "population_size": population_size,
            "chosen_families": len(chosen_fam_ids),
            "realized_vaf": realized_vaf,
            "complete_duplex": complete_duplex_cnt,
            "single_strand_duplex": single_strand_duplex_cnt,
            "unparseable_tags": unparseable_tags_count,
            "status": "ACTIVE",
            "selection_rule": selection_rule,
            "bg_total_families": len(bg_fam_counts),
            "bg_alt_families": sum(1 for _t, _a in bg_fam_counts.values() if _a > 0 and _a * 2 >= _t),
            **sb_info,
            **ob_info,
            **cons_info,
            **comp_info
        }

    bam.close()
    return locus_reads, read_boost_map, variant_stats


_WARNED_NON_SNV = set()


def mutate_alignment_record_multi(r_obj, active_mutations, umi_opts, bg_model, aligned_pairs=None):
    """Applies SNV substitutions to an alignment record.

    Returns (record, mutated_indices, injected_by_variant, ref_mismatch_ids):
      injected_by_variant: {cosmic_id: [injected_qual, ...]} for variants that were ACTUALLY applied.
      ref_mismatch_ids:    set of cosmic_ids skipped because the read base did not equal the manifest ref.
    Non-SNVs and variants whose position is not aligned in this read appear in neither output.
    """
    if not active_mutations:
        return r_obj, set(), {}, set()

    seq_list = list(r_obj.query_sequence)
    qual_list = list(r_obj.query_qualities) if r_obj.query_qualities is not None else [40] * len(seq_list)

    if len(qual_list) < len(seq_list):
        qual_list = qual_list + [40] * (len(seq_list) - len(qual_list))
    elif len(qual_list) > len(seq_list):
        qual_list = qual_list[:len(seq_list)]

    sorted_muts = sorted(active_mutations, key=lambda x: x["pos"], reverse=True)
    mutated_indices = set()
    injected_by_variant = {}
    ref_mismatch_ids = set()

    if aligned_pairs is None:
        aligned_pairs = r_obj.get_aligned_pairs(matches_only=True) if not r_obj.is_unmapped else []

    qual_mode = umi_opts.get("qual_mode", "keep")
    manual_qual = umi_opts.get("manual_qual", 30)

    for mut in sorted_muts:
        vid = mut["cosmic_id"]
        ref_base = mut.get("ref", "").strip().upper()
        alt_base = mut["alt"].strip().upper()

        if len(ref_base) != 1 or len(alt_base) != 1:
            if vid not in _WARNED_NON_SNV:  # warn once per variant, not once per read
                _WARNED_NON_SNV.add(vid)
                print(f"[WARNING] Non-SNV variant encountered ({vid}: {ref_base}->{alt_base}). Skipping mutation.")
            continue

        target_ref_pos = mut["pos"] - 1

        for read_pos, ref_pos in aligned_pairs:
            if ref_pos == target_ref_pos and read_pos is not None and read_pos < len(seq_list):
                # query_sequence is stored in reference orientation in BAM, and aligned_pairs
                # excludes soft-clipped bases, so this comparison is orientation-safe.
                if seq_list[read_pos].upper() != ref_base:
                    ref_mismatch_ids.add(vid)
                    break

                seq_list[read_pos] = alt_base
                mutated_indices.add(read_pos)

                cycle_pos = (len(seq_list) - 1 - read_pos) if r_obj.is_reverse else read_pos

                if qual_mode == "decay" and bg_model is not None:
                    q_char = calculate_injected_qual_v2(chr(qual_list[read_pos] + 33), cycle_pos, bg_model)
                    final_q = ord(q_char) - 33
                elif qual_mode == "manual":
                    final_q = manual_qual
                else:  # keep mode
                    final_q = qual_list[read_pos]

                qual_list[read_pos] = final_q
                injected_by_variant[vid] = [final_q]
                break

    if not mutated_indices:
        return r_obj, set(), {}, ref_mismatch_ids

    r_obj.query_sequence = "".join(seq_list)
    r_obj.query_qualities = qual_list
    return r_obj, mutated_indices, injected_by_variant, ref_mismatch_ids


def _quals_to_phred_bytes(query_qualities, seq_len):
    """Direct byte array conversion for zero-copy high throughput processing."""
    if query_qualities is None:
        return b"I" * seq_len

    q_bytes = bytes(query_qualities)
    if len(q_bytes) < seq_len:
        q_bytes = q_bytes + b"\x28" * (seq_len - len(q_bytes))
    elif len(q_bytes) > seq_len:
        q_bytes = q_bytes[:seq_len]

    return q_bytes.translate(PHRED_TRANS_TABLE)


def inject_noise_to_record(r_obj, norm_chrom, bg_model, error_opts, mutated_indices=None, aligned_pairs=None, py_rng=None):
    seq = r_obj.query_sequence
    quals = list(r_obj.query_qualities) if r_obj.query_qualities is not None else [40] * len(seq)
    qual_str = _quals_to_phred_bytes(quals, len(seq)).decode("ascii")

    scale = error_opts.get("scale", 1.0)
    flat_rate_v2 = error_opts.get("flat_rate") if error_opts.get("flat_rate") is not None else 0.0005
    flat_rate_v1 = error_opts.get("flat_rate") if error_opts.get("flat_rate") is not None else 0.001

    # Fallback aligned_pairs to None if empty or unmapped to allow r_obj alignment recovery in apply_sequencing_error_model_v2
    pass_aligned_pairs = aligned_pairs if (aligned_pairs and len(aligned_pairs) > 0) else None

    if bg_model is not None and hasattr(bg_model, "position_probs") and bg_model.position_probs:
        seq = apply_sequencing_error_model_v2(
            seq, qual_str, norm_chrom, r_obj.reference_start, bg_model,
            mutated_indices, scale, flat_rate_v2, r_obj=r_obj, aligned_pairs=pass_aligned_pairs, py_rng=py_rng
        )
    else:
        seq = apply_sequencing_error_model(seq, qual_str, mutated_indices, scale, flat_rate_v1, py_rng=py_rng)

    r_obj.query_sequence = seq
    if len(quals) != len(seq):
        quals = quals[:len(seq)] if len(quals) > len(seq) else quals + [40] * (len(seq) - len(quals))
    r_obj.query_qualities = quals
    return r_obj


def generate_decaying_qual_string(length, start_pos=0, bg_model=None, base_qual=38, min_qual=20):
    if length <= 0:
        return ""

    cycles = np.arange(start_pos, start_pos + length, dtype=np.int32)

    if bg_model is not None and hasattr(bg_model, "cycle_qual_array") and len(bg_model.cycle_qual_array) > 0:
        valid_mask = cycles < len(bg_model.cycle_qual_array)
        valid_cycles = cycles[valid_mask]
        means = np.full(length, base_qual, dtype=np.float32)

        if np.any(valid_mask):
            lookup_vals = bg_model.cycle_qual_array[valid_cycles]
            found_mask = lookup_vals >= 0
            means[np.flatnonzero(valid_mask)[found_mask]] = lookup_vals[found_mask]

        q_vals = np.clip(np.random.normal(loc=means, scale=2.5), min_qual, 41)
    else:
        decay_factors = np.exp(-0.015 * cycles)
        q_vals = np.maximum(min_qual, np.round(base_qual * decay_factors))

    q_bytes = (q_vals.astype(np.uint8) + 33).tobytes()
    return q_bytes.decode("ascii")


def format_read_as_fastq_bytes(read, read_num, umi_opts, duplex_opts, bg_model=None, tracker=None, static_umi_qual_b=b"", static_spacer_qual_b=b""):
    seq_str = read.query_sequence
    qual_bytes = _quals_to_phred_bytes(read.query_qualities, len(seq_str))

    if read.is_reverse:
        seq_str = reverse_complement(seq_str)
        qual_bytes = qual_bytes[::-1]

    seq_b = seq_str.encode("ascii")
    qual_b = qual_bytes

    placement = umi_opts.get("umi_placement", "inline")
    spacer_len = umi_opts.get("spacer_len") or 0
    umi_len = umi_opts.get("umi_len", 12)
    fallback_umi_b = b"N" * umi_len
    kit_spacer_b = b"N" * spacer_len

    attach_umi = (duplex_opts.get("mode") == "duplex") or (read_num == umi_opts.get("umi_read"))

    if attach_umi:
        raw_umi = extract_umi_sequence(read, read_num, umi_len)
        if not raw_umi:
            raw_umi_b = fallback_umi_b
            is_fallback = True
        else:
            raw_umi_b = raw_umi.encode("ascii")
            is_fallback = False

        if tracker:
            tracker.log_result(is_fallback)
    else:
        raw_umi_b = b""
        is_fallback = False

    if attach_umi and placement == "inline":
        if static_umi_qual_b and static_spacer_qual_b:
            umi_qual_b = static_umi_qual_b
            spacer_qual_b = static_spacer_qual_b
        else:
            umi_qual_b = generate_decaying_qual_string(len(raw_umi_b), 0, bg_model).encode("ascii")
            spacer_qual_b = generate_decaying_qual_string(spacer_len, len(raw_umi_b), bg_model).encode("ascii")

        seq_b = b"".join([raw_umi_b, kit_spacer_b, seq_b])
        qual_b = b"".join([umi_qual_b, spacer_qual_b, qual_b])

    qname_b = sanitize_qname(read.query_name).encode("ascii")
    mate_comment_b = b" 1:N:0:1\n" if read_num == 1 else b" 2:N:0:1\n"

    if placement == "header":
        # Same token on BOTH mates so R1/R2 read names stay identical.
        token = resolve_header_umi_token(read, umi_opts, duplex_opts) or ("N" * umi_len)
        header_b = b"@" + qname_b + b":" + token.encode("ascii") + mate_comment_b
    else:
        header_b = b"@" + qname_b + mate_comment_b

    return b"".join([header_b, seq_b, b"\n+\n", qual_b, b"\n"])


def write_spike_ledger(ledger_path, active_variants, ledger_counts, variant_stats):
    """Writes detailed CSV spike ledger for auditing exact variant injection metrics."""
    fieldnames = [
        "cosmic_id", "chrom", "pos", "ref", "alt", "target_vaf", "selection_rule",
        "eligible_families", "chosen_families", "realized_family_vaf_pct",
        "overlapping_reads_r1", "overlapping_reads_r2",
        "strand_bias_skipped_reads", "ref_mismatch_skipped", "reads_actually_mutated", "boosted_copies_written",
        "orphans_emitted", "mutant_reads_written", "total_reads_written",
        "expected_consensus_vaf_pct", "consensus_units_total", "consensus_units_alt",
        "consensus_units_discordant", "consensus_alt_fwd_fraction", "avg_injected_base_qual",
        "requested_strand_bias", "sb_reachable_min", "sb_reachable_max",
        "alt_reads_fwd_written", "alt_reads_rev_written", "realized_fwd_alt_fraction",
        "strand_bias_status", "sb_alt_units_potential", "sb_alt_units_planned",
        "pool_fwd_only_families", "pool_rev_only_families", "pool_both_families",
        "background_alt_read_pct", "background_nonref_read_pct", "background_alt_family_pct",
        "requested_orientation_bias", "ob_reachable_min", "ob_reachable_max", "alt_families_strand_A",
        "alt_families_strand_B", "realized_strand_A_fraction", "orientation_bias_status",
        "ob_forced_A_families", "ob_forced_B_families", "ob_free_families",
        "chosen_fwd_only_families", "chosen_rev_only_families", "chosen_both_families",
        "pool_F1R2_families", "pool_F2R1_families", "pool_both_pairtype_families",
        "chosen_F1R2_families", "chosen_F2R1_families", "chosen_both_pairtype_families",
        "native_F1R2_read_pct", "chosen_F1R2_read_pct", "alt_reads_F1R2_written", "alt_reads_F2R1_written",
        "chosen_r1_5p_dist_min", "chosen_r1_5p_dist_median", "chosen_r2_5p_dist_min", "chosen_r2_5p_dist_median",
        "native_r1_5p_dist_min", "native_r1_5p_dist_median", "native_r2_5p_dist_min", "native_r2_5p_dist_median",
        "bed_region", "min_edge_dist_bp"
    ]

    with open(ledger_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()

        for v in active_variants:
            vid = v["cosmic_id"]
            st = variant_stats[vid]
            counts = ledger_counts[vid]

            mut_reads = counts["reads_actually_mutated"]
            boosted_writes = counts["boosted_copies_written"]
            mut_written = counts["mutant_reads_written"]
            tot_written = counts["total_reads_written"]

            # (mutated reads x copies) / (all overlapping reads x copies); strand-bias-skipped and
            # ref-mismatch reads stay in the denominator.
            avg_qual = (sum(counts["injected_quals"]) / len(counts["injected_quals"])) if counts["injected_quals"] else 0.0

            def _fmt(x):
                return "NA" if x is None else f"{x:.4f}"

            def _bg_pct(num, den):
                # background of the unspiked reads at the site; NA when it could not be measured
                return "NA" if (num is None or not den) else f"{num / den * 100:.3f}%"

            def _pt_pct(n1, n2):
                # share of reads covering the site that are F1R2 (read 1 forward); NA when no paired reads were seen
                return "NA" if (n1 is None or n2 is None or (n1 + n2) == 0) else f"{n1 / (n1 + n2) * 100:.1f}%"

            def _ds_fmt(d, rn, k):
                if not d or d.get(rn) is None:
                    return "NA"
                x = d[rn][k]
                return str(int(x)) if float(x).is_integer() else f"{x:.1f}"

            n_fwd, n_rev = counts["alt_reads_fwd_written"], counts["alt_reads_rev_written"]
            # status is judged on consensus units (what the caller sees); read-level columns stay for reference
            realized_fwd = st.get("cons_alt_fwd_fraction")
            n_units = st.get("cons_units_alt", 0)
            if not st.get("sb_enabled"):
                sb_status = "NA"
            elif st.get("sb_plan_status") == "UNREACHABLE":
                sb_status = "UNREACHABLE"
            elif realized_fwd is None:
                sb_status = "NO_ALT_READS"
            elif abs(realized_fwd - st["sb_requested"]) <= max(0.02, 1.0 / max(1, n_units)):
                sb_status = "MET"
            else:
                sb_status = "APPROX"

            writer.writerow({
                "cosmic_id": vid,
                "chrom": v["chrom"],
                "pos": v["pos"],
                "ref": v["ref"],
                "alt": v["alt"],
                "target_vaf": f"{v['vaf']*100:.3f}%",
                "selection_rule": st.get("selection_rule", "binomial"),
                "eligible_families": st["population_size"],
                "chosen_families": st["chosen_families"],
                "realized_family_vaf_pct": f"{st['realized_vaf']:.3f}%",
                "overlapping_reads_r1": counts["overlapping_reads_r1"],
                "overlapping_reads_r2": counts["overlapping_reads_r2"],
                "strand_bias_skipped_reads": counts["strand_bias_skipped_reads"],
                "ref_mismatch_skipped": counts["ref_mismatch_skipped"],
                "reads_actually_mutated": mut_reads,
                "boosted_copies_written": boosted_writes,
                "orphans_emitted": counts["orphans_emitted"],
                "mutant_reads_written": mut_written,
                "total_reads_written": tot_written,
                "expected_consensus_vaf_pct": f"{st.get('expected_consensus_vaf', 0.0):.3f}%",
                "consensus_units_total": st.get("cons_units_total", 0),
                "consensus_units_alt": st.get("cons_units_alt", 0),
                "consensus_units_discordant": st.get("cons_units_discordant", 0),
                "consensus_alt_fwd_fraction": ("NA" if st.get("cons_alt_fwd_fraction") is None else f"{st['cons_alt_fwd_fraction']:.4f}"),
                "avg_injected_base_qual": f"{avg_qual:.1f}",
                "requested_strand_bias": _fmt(st["sb_requested"]) if st.get("sb_enabled") else "NA",
                "sb_reachable_min": _fmt(st.get("sb_min")) if st.get("sb_enabled") else "NA",
                "sb_reachable_max": _fmt(st.get("sb_max")) if st.get("sb_enabled") else "NA",
                "alt_reads_fwd_written": n_fwd,
                "alt_reads_rev_written": n_rev,
                "realized_fwd_alt_fraction": _fmt((n_fwd / (n_fwd + n_rev)) if (n_fwd + n_rev) > 0 else None),
                "strand_bias_status": sb_status,
                "sb_alt_units_potential": st.get("sb_alt_units_potential", "NA"),
                "sb_alt_units_planned": st.get("sb_alt_units_planned", "NA"),
                "pool_fwd_only_families": st.get("sb_pool_fwd_only", "NA"),
                "pool_rev_only_families": st.get("sb_pool_rev_only", "NA"),
                "pool_both_families": st.get("sb_pool_both", "NA"),
                "background_alt_read_pct": _bg_pct(v.get("bg_reads", {}).get("alt_reads"), v.get("bg_reads", {}).get("total_reads")),
                "background_nonref_read_pct": _bg_pct(v.get("bg_reads", {}).get("nonref_reads"), v.get("bg_reads", {}).get("total_reads")),
                "background_alt_family_pct": _bg_pct(st.get("bg_alt_families"), st.get("bg_total_families")),
                "requested_orientation_bias": _fmt(st["ob_requested"]) if st.get("ob_enabled") else "NA",
                "ob_reachable_min": _fmt(st.get("ob_min")) if st.get("ob_enabled") else "NA",
                "ob_reachable_max": _fmt(st.get("ob_max")) if st.get("ob_enabled") else "NA",
                "alt_families_strand_A": st.get("ob_families_a", "NA"),
                "alt_families_strand_B": st.get("ob_families_b", "NA"),
                "realized_strand_A_fraction": _fmt(st.get("ob_planned_a")) if st.get("ob_enabled") else "NA",
                "orientation_bias_status": st.get("ob_status", "NA") if st.get("ob_enabled") else "NA",
                "ob_forced_A_families": st.get("ob_forced_a", "NA"),
                "ob_forced_B_families": st.get("ob_forced_b", "NA"),
                "ob_free_families": st.get("ob_free", "NA"),
                "chosen_fwd_only_families": st.get("chosen_fwd_only", "NA"),
                "chosen_rev_only_families": st.get("chosen_rev_only", "NA"),
                "chosen_both_families": st.get("chosen_both", "NA"),
                "pool_F1R2_families": st.get("pool_f1r2_fams", "NA"),
                "pool_F2R1_families": st.get("pool_f2r1_fams", "NA"),
                "pool_both_pairtype_families": st.get("pool_both_pt_fams", "NA"),
                "chosen_F1R2_families": st.get("chosen_f1r2_fams", "NA"),
                "chosen_F2R1_families": st.get("chosen_f2r1_fams", "NA"),
                "chosen_both_pairtype_families": st.get("chosen_both_pt_fams", "NA"),
                "native_F1R2_read_pct": _pt_pct(st.get("native_f1r2_reads"), st.get("native_f2r1_reads")),
                "chosen_F1R2_read_pct": _pt_pct(st.get("chosen_f1r2_reads"), st.get("chosen_f2r1_reads")),
                "alt_reads_F1R2_written": counts["alt_reads_f1r2_written"],
                "alt_reads_F2R1_written": counts["alt_reads_f2r1_written"],
                "chosen_r1_5p_dist_min": _ds_fmt(st.get("chosen_dist"), 1, 0),
                "chosen_r1_5p_dist_median": _ds_fmt(st.get("chosen_dist"), 1, 1),
                "chosen_r2_5p_dist_min": _ds_fmt(st.get("chosen_dist"), 2, 0),
                "chosen_r2_5p_dist_median": _ds_fmt(st.get("chosen_dist"), 2, 1),
                "native_r1_5p_dist_min": _ds_fmt(st.get("native_dist"), 1, 0),
                "native_r1_5p_dist_median": _ds_fmt(st.get("native_dist"), 1, 1),
                "native_r2_5p_dist_min": _ds_fmt(st.get("native_dist"), 2, 0),
                "native_r2_5p_dist_median": _ds_fmt(st.get("native_dist"), 2, 1),
                "bed_region": v.get("bed_region_tag", "NA") if v.get("bed_annotated") else "NA",
                "min_edge_dist_bp": (v["min_edge_dist"] if (v.get("bed_annotated") and v.get("min_edge_dist", -1) >= 0) else "NA")
            })

    print(f"[+] Spike ledger report written to '{ledger_path}'.")


def run_global_bam_to_fastq_multi(
    coord_bam, out_r1_path, out_r2_path, variants_list, umi_opts, boost_opts, error_opts, duplex_opts,
    orphan_opts=None, ledger_path=None, mutated_reads_path=None, bg_model=None, bam_threads=4, global_seed=42,
    strict_strand_bias=False,
    strict_orientation_bias=False,
    no_spike=False
):
    if orphan_opts is None:
        orphan_opts = {"mode": "separate", "out_orphan_r1": f"{out_r1_path}.orphans.fastq.gz", "out_orphan_r2": f"{out_r2_path}.orphans.fastq.gz"}

    validate_bam_sort_and_index(coord_bam)

    variants_by_id = {v["cosmic_id"]: v for v in variants_list}

    noise_threshold = umi_opts.get("noise_threshold", 0.01)
    bg_warn_fraction = umi_opts.get("background_warn_fraction", 0.2)
    for v in variants_list:
        bg = {}
        depth = check_bam_integrity_and_variants(coord_bam, v["chrom"], v["pos"], v["ref"],
                                                 force_noise=umi_opts.get("force_noise", False),
                                                 alt_base=v.get("alt"), noise_threshold=noise_threshold, bg=bg)
        v["depth"] = depth
        v["bg_reads"] = bg
        if depth == -1:
            v["excluded_due_to_noise"] = True
        # Relative warning (never excludes): background reads that already carry the spiked alt add to the final VAF.
        if bg_warn_fraction > 0 and bg.get("total_reads", 0) > 0 and v.get("vaf", 0) > 0:
            bg_alt_frac = bg["alt_reads"] / bg["total_reads"]
            if bg_alt_frac > bg_warn_fraction * v["vaf"]:
                print(f"[WARNING] {v['cosmic_id']} ({v['chrom']}:{v['pos']}): {bg['alt_reads']:,} of {bg['total_reads']:,} unspiked reads "
                      f"({bg_alt_frac*100:.3f}%) already carry the spiked alt, which is {bg_alt_frac / v['vaf'] * 100:.0f}% of the target VAF "
                      f"({v['vaf']*100:.3f}%). The observed VAF will be about target plus background.")

    if no_spike:
        if variants_list:
            sys.exit("[FATAL] Internal error: --no-spike with a non-empty variant list.")
        # Nothing to scan: every read takes the untouched-read branch of pass 2.
        locus_reads, read_boost_map, variant_stats = {}, {}, {}
    else:
        locus_reads, read_boost_map, variant_stats = scan_target_families_multi(
            coord_bam, variants_list, boost_opts=boost_opts, duplex_opts=duplex_opts,
            stochastic=umi_opts.get("stochastic", True), global_seed=global_seed
        )

    if not no_spike:
        print("\nVariant Pre-Flight Diagnostics:")
        print("-" * 109)
        print(f"{'COSMIC ID':<14} {'Locus':<15} {'Depth':<8} {'5pMin(sel)':<11} {'5pMed(all)':<11} {'Comp-Duplex':<13} {'Single-Duplex':<15} {'Unparseable':<12} {'Target VAF':<10}")
        print("-" * 109)

        for v in variants_list:
            vid = v["cosmic_id"]
            st = variant_stats[vid]
            locus_str = f"{v['chrom']}:{v['pos']}"
            depth_str = f"{v.get('depth', 0):,}" if v.get('depth', 0) >= 0 else "EXCL"
            # Distance from the site to the read's 5' end (aligned bases): smallest over the selected families' reads
            # (both reads) and median over all unspiked reads at the site. Read-based, not the BED edge.
            def _rd(d, k):
                if not d or d.get("all") is None:
                    return "N/A"
                x = d["all"][k]
                return f"{int(x)} bp" if float(x).is_integer() else f"{x:.1f} bp"
            edge_min_str = _rd(st.get("chosen_dist"), 0)
            edge_med_str = _rd(st.get("native_dist"), 1)
            comp_dup = f"{st['complete_duplex']:,}"
            sing_dup = f"{st['single_strand_duplex']:,}"
            unparseable = f"{st['unparseable_tags']:,}"
            target_vaf_str = f"{v['vaf']*100:.2f}%"

            print(f"{vid:<14} {locus_str:<15} {depth_str:<8} {edge_min_str:<11} {edge_med_str:<11} {comp_dup:<13} {sing_dup:<15} {unparseable:<12} {target_vaf_str:<10}")
        print("-" * 109 + "\n")

    active_variants = [v for v in variants_list if variant_stats[v["cosmic_id"]]["status"] == "ACTIVE"]

    if not active_variants and not no_spike:
        sys.exit("[FATAL] All manifest variants were excluded due to zero family coverage, BED filtering, or background noise.")

    sb_variants = [v for v in active_variants if variant_stats[v["cosmic_id"]].get("sb_enabled")]
    if sb_variants:
        print("Strand-Bias Plan (fraction of alt reads that are forward-aligned):")
        print("---------------------------------------------------------------------------------------------")
        print(f"{'COSMIC ID':<14} {'Requested':<10} {'Reachable range':<17} {'Planned':<9} {'Status':<13} {'Alt units kept':<15} {'Pool F/R/both'}")
        print("---------------------------------------------------------------------------------------------")
        unreachable = []
        for v in sb_variants:
            st = variant_stats[v["cosmic_id"]]
            lo, hi, pl = st.get("sb_min"), st.get("sb_max"), st.get("sb_planned_fwd")
            rng_s = f"{lo:.2f}-{hi:.2f}" if (lo is not None and hi is not None) else "NA"
            pl_s = f"{pl:.3f}" if pl is not None else "NA"
            pool_s = f"{st['sb_pool_fwd_only']}/{st['sb_pool_rev_only']}/{st['sb_pool_both']}"
            print(f"{v['cosmic_id']:<14} {st['sb_requested']:<10.2f} {rng_s:<17} {pl_s:<9} {st['sb_plan_status']:<13} {str(st['sb_alt_units_planned']) + '/' + str(st['sb_alt_units_potential']):<15} {pool_s}")
            if st["sb_plan_status"] == "UNREACHABLE":
                unreachable.append(v["cosmic_id"])
                nearest = lo if st["sb_requested"] < lo else hi
                print(f"  [WARNING] {v['cosmic_id']}: requested strand bias {st['sb_requested']:.2f} is outside what the "
                      f"covering reads allow ({rng_s}); the closest reachable value (~{nearest:.2f}) will be used. "
                      f"Most chosen families have a consensus read in one orientation only at this site.")
        print("---------------------------------------------------------------------------------------------\n")
        if strict_strand_bias and unreachable:
            sys.exit(f"[FATAL] --strict-strand-bias: unreachable strand-bias request for {', '.join(unreachable)}.")

    ob_variants = [v for v in active_variants if variant_stats[v["cosmic_id"]].get("ob_enabled")]
    if ob_variants:
        print("Orientation-Bias Plan (fraction of alt-carrying families with the alt on strand A, /A = F1R2 in duplex data):")
        print("---------------------------------------------------------------------------------------------")
        print(f"{'COSMIC ID':<14} {'Requested':<10} {'Reachable range':<17} {'Planned':<9} {'Status':<13} {'Families A/B':<14} {'Forced A/B, free'}")
        print("---------------------------------------------------------------------------------------------")
        ob_unreachable = []
        for v in ob_variants:
            st = variant_stats[v["cosmic_id"]]
            lo, hi, pl = st.get("ob_min"), st.get("ob_max"), st.get("ob_planned_a")
            rng_s = f"{lo:.2f}-{hi:.2f}" if (lo is not None and hi is not None) else "NA"
            pl_s = f"{pl:.3f}" if pl is not None else "NA"
            fam_s = f"{st['ob_families_a']}/{st['ob_families_b']}"
            forced_s = f"{st['ob_forced_a']}/{st['ob_forced_b']}, {st['ob_free']}"
            print(f"{v['cosmic_id']:<14} {st['ob_requested']:<10.2f} {rng_s:<17} {pl_s:<9} {st['ob_status']:<13} {fam_s:<14} {forced_s}")
            if st["ob_status"] == "UNREACHABLE":
                ob_unreachable.append(v["cosmic_id"])
                nearest = lo if st["ob_requested"] < lo else hi
                print(f"  [WARNING] {v['cosmic_id']}: requested orientation bias {st['ob_requested']:.2f} is outside what the chosen "
                      f"families allow ({rng_s}); the closest reachable value (~{nearest:.2f}) will be used. "
                      f"{st['ob_forced_a']} chosen families are covered by strand A only and {st['ob_forced_b']} by strand B only at this site, "
                      f"so their strand cannot be changed.")
            elif st["ob_status"] == "APPROX":
                print(f"  [WARNING] {v['cosmic_id']}: orientation bias planned at {pl_s} instead of the requested "
                      f"{st['ob_requested']:.2f} (too few families to hit it exactly).")
            elif st["ob_status"] == "NO_ALT_FAMILIES":
                print(f"  [WARNING] {v['cosmic_id']}: no alt-carrying families, orientation bias not applied.")
        print("---------------------------------------------------------------------------------------------\n")
        if strict_orientation_bias and ob_unreachable:
            sys.exit(f"[FATAL] --strict-orientation-bias: unreachable orientation-bias request for {', '.join(ob_unreachable)}.")

    if no_spike:
        print("\n" + "="*70 + "\n NO-SPIKE MODE: BASELINE FASTQs, NO VARIANTS SPIKED\n" + "="*70)
    else:
        print(f"\n" + "="*70 + f"\n MULTI-SNV DIGITAL SPIKE-IN TRANSFORMATION SUMMARY\n" + "="*70)
    print(f" Protocol Mode:                 {duplex_opts['mode'].upper()}")
    print(f" UMI Placement Mode:             {umi_opts.get('umi_placement', 'inline').upper()}")
    if duplex_opts["mode"] == "duplex":
        print(f" Duplex Strand Source:           {duplex_opts['strand_source']}")
    print(f" Target Variants Count:          {len(variants_list)} ({len(active_variants)} Active)")
    if no_spike:
        print(" Reads are written unmodified; no ledger or mutated-read file is produced.")
    print("-" * 70)
    for v in variants_list:
        vid = v["cosmic_id"]
        st = variant_stats[vid]
        if st["status"] == "EXCLUDED_ZERO_COVERAGE":
            print(f" Variant {vid:15s} | Status: EXCLUDED (0 Coverage)")
        elif st["status"] == "EXCLUDED_NOISE":
            print(f" Variant {vid:15s} | Status: EXCLUDED (Background Noise > 1%)")
        else:
            print(f" Variant {vid:15s} | Target VAF: {v['vaf']*100:6.3f}% | Expected consensus VAF: {st['expected_consensus_vaf']:6.3f}% | Families chosen: {st['realized_vaf']:6.3f}% of {st['population_size']:,}")
    print("="*70 + "\n")

    total_cpu_cores = os.cpu_count() or 4
    budget = max(1, total_cpu_cores // 3)
    reader_threads = min(bam_threads, budget)
    pigz_threads = min(bam_threads, budget)
    
    print(f"[-] Initializing FASTQ Writers (R1/R2 pigz threads: {pigz_threads} each, BAM reader threads: {reader_threads})...")
    print(f"[-] Pass 2: Streaming full BAM sequence data directly to paired FASTQs...")

    with ExitStack() as stack:
        out_r1 = stack.enter_context(FastqChunkWriter(out_r1_path, threads=pigz_threads))
        out_r2 = stack.enter_context(FastqChunkWriter(out_r2_path, threads=pigz_threads))

        orphan_mode = orphan_opts.get("mode", "separate")
        if orphan_mode == "separate":
            out_orph_r1 = stack.enter_context(FastqChunkWriter(orphan_opts["out_orphan_r1"], threads=pigz_threads))
            out_orph_r2 = stack.enter_context(FastqChunkWriter(orphan_opts["out_orphan_r2"], threads=pigz_threads))
        else:
            out_orph_r1, out_orph_r2 = None, None

        umi_tracker = UmiEmbedTracker()
        bam = stack.enter_context(pysam.AlignmentFile(coord_bam, "rb", threads=reader_threads))
        bam_refs = set(bam.references)

        err_enabled = error_opts.get("enabled", False)
        force_global_err = error_opts.get("force_global", False)
        boost_enabled = boost_opts.get("enabled", False)

        placement = umi_opts.get("umi_placement", "inline")
        spacer_len = umi_opts.get("spacer_len") or 0
        umi_len = umi_opts.get("umi_len", 12)
        fallback_umi = "N" * umi_len
        kit_spacer = "N" * spacer_len

        has_empirical_qual = bg_model is not None and hasattr(bg_model, "cycle_qual_array") and len(bg_model.cycle_qual_array) > 0
        if not has_empirical_qual:
            static_umi_qual_b = generate_decaying_qual_string(umi_len, 0, bg_model=None).encode("ascii")
            static_spacer_qual_b = generate_decaying_qual_string(spacer_len, umi_len, bg_model=None).encode("ascii")
        else:
            static_umi_qual_b = b""
            static_spacer_qual_b = b""

        duplex_mode = duplex_opts.get("mode")
        umi_read_target = umi_opts.get("umi_read")

        # Metrics for spike ledger
        ledger_counts = defaultdict(lambda: {
            "overlapping_reads_r1": 0, "overlapping_reads_r2": 0,
            "strand_bias_skipped_reads": 0, "ref_mismatch_skipped": 0, "reads_actually_mutated": 0,
            "boosted_copies_written": 0, "orphans_emitted": 0, "injected_quals": [],
            "mutant_reads_written": 0, "total_reads_written": 0,
            "alt_reads_fwd_written": 0, "alt_reads_rev_written": 0,
            "alt_reads_f1r2_written": 0, "alt_reads_f2r1_written": 0
        })
        mutated_qnames_set = set()

        def process_and_emit_read(read_obj, is_orphan=False, target_writer=None, forced_boost_count=None):
            qname = read_obj.query_name
            read_num = 1 if read_obj.is_read1 else 2

            if target_writer is not None:
                out_writer = target_writer
            elif is_orphan:
                if orphan_mode == "discard":
                    return
                out_writer = out_orph_r1 if read_num == 1 else out_orph_r2
            else:
                out_writer = out_r1 if read_num == 1 else out_r2

            lookup_qname = strip_boost_suffix(qname)

            if lookup_qname not in locus_reads:
                rec_bytes = format_read_as_fastq_bytes(
                    read_obj, read_num, umi_opts, duplex_opts, bg_model, 
                    tracker=umi_tracker, static_umi_qual_b=static_umi_qual_b, static_spacer_qual_b=static_spacer_qual_b
                )
                out_writer.write_record(rec_bytes)
                return

            read_var_map = locus_reads[lookup_qname]
            is_any_mutant = any(v["is_mutant"] for v in read_var_map.values())
            aligned_pairs = read_obj.get_aligned_pairs(matches_only=True) if (not read_obj.is_unmapped) else []

            repeat_count = forced_boost_count if forced_boost_count is not None else (read_boost_map.get(lookup_qname, 1) if boost_enabled else 1)

            active_mutations = []
            seq = read_obj.query_sequence

            if seq and not read_obj.is_unmapped:
                for vid, var_status in read_var_map.items():
                    v_info = variants_by_id[vid]
                    target_ref_pos = v_info["pos"] - 1

                    overlaps_pos = any(ref_p == target_ref_pos and read_p is not None and read_p < len(seq)
                                       for read_p, ref_p in aligned_pairs)

                    if overlaps_pos:
                        if read_num == 1:
                            ledger_counts[vid]["overlapping_reads_r1"] += 1
                        else:
                            ledger_counts[vid]["overlapping_reads_r2"] += 1
                        ledger_counts[vid]["total_reads_written"] += repeat_count

                        if var_status["is_mutant"]:
                            # Strand bias is planned per family in pass 1: 'mut_orient' is the aligned orientation
                            # ('F'/'R') that carries the alt, or None when every covering read is mutated.
                            mut_orient = var_status.get("mut_orient")
                            if mut_orient is None or mut_orient == read_direction(read_obj.is_reverse):
                                active_mutations.append(v_info)
                            else:
                                ledger_counts[vid]["strand_bias_skipped_reads"] += 1

            # Per-locus py_rng for noise injection
            var_noise_seed = int.from_bytes(hashlib.sha256(f"{global_seed}:{qname}:{read_num}".encode("utf-8")).digest()[:4], byteorder="big")
            noise_py_rng = random.Random(var_noise_seed)

            r_clean, mut_indices, injected_by_variant, ref_mismatch_ids = mutate_alignment_record_multi(
                read_obj, active_mutations, umi_opts, bg_model, aligned_pairs=aligned_pairs
            )

            for vid in ref_mismatch_ids:
                ledger_counts[vid]["ref_mismatch_skipped"] += 1

            if mut_indices:
                mutated_qnames_set.add(sanitize_qname(qname))
                # Only variants that were actually applied to this read are credited.
                for vid, v_quals in injected_by_variant.items():
                    c = ledger_counts[vid]
                    c["reads_actually_mutated"] += 1
                    c["mutant_reads_written"] += repeat_count      # original + its boosted copies
                    c["alt_reads_rev_written" if read_obj.is_reverse else "alt_reads_fwd_written"] += repeat_count
                    if read_obj.is_read1 or read_obj.is_read2:
                        c["alt_reads_f1r2_written" if (read_obj.is_read1 != read_obj.is_reverse) else "alt_reads_f2r1_written"] += repeat_count
                    c["boosted_copies_written"] += (repeat_count - 1)
                    c["injected_quals"].extend(v_quals)
                    if is_orphan:
                        c["orphans_emitted"] += repeat_count

            apply_orig_noise = err_enabled and (bool(mut_indices) or force_global_err)
            apply_dup_noise = err_enabled and (repeat_count > 1)

            base_seq = str(r_clean.query_sequence)
            base_qual = list(r_clean.query_qualities) if r_clean.query_qualities is not None else [40] * len(base_seq)

            ref_chrom = normalize_chrom(r_clean.reference_name, bam_refs) if r_clean.reference_name else None

            if apply_orig_noise:
                orig_read = inject_noise_to_record(r_clean, ref_chrom, bg_model, error_opts, mut_indices, aligned_pairs=aligned_pairs, py_rng=noise_py_rng)
            else:
                orig_read = r_clean

            out_writer.write_record(format_read_as_fastq_bytes(
                orig_read, read_num, umi_opts, duplex_opts, bg_model, 
                tracker=umi_tracker, static_umi_qual_b=static_umi_qual_b, static_spacer_qual_b=static_spacer_qual_b
            ))

            if repeat_count > 1:
                raw_umi = extract_umi_sequence(orig_read, read_num, umi_len)
                raw_umi_b = raw_umi.encode("ascii") if raw_umi else fallback_umi.encode("ascii")
                attach_umi = (duplex_mode == "duplex") or (read_num == umi_read_target)

                if attach_umi and placement == "inline":
                    if static_umi_qual_b and static_spacer_qual_b:
                        umi_qual_b = static_umi_qual_b
                        spacer_qual_b = static_spacer_qual_b
                    else:
                        umi_qual_b = generate_decaying_qual_string(len(raw_umi_b), 0, bg_model).encode("ascii")
                        spacer_qual_b = generate_decaying_qual_string(spacer_len, len(raw_umi_b), bg_model).encode("ascii")

                    kit_spacer_b = kit_spacer.encode("ascii")
                    boost_seq_prefix_b = b"".join([raw_umi_b, kit_spacer_b])
                    boost_qual_prefix_b = b"".join([umi_qual_b, spacer_qual_b])
                else:
                    boost_seq_prefix_b = b""
                    boost_qual_prefix_b = b""

                final_qname = sanitize_qname(qname)
                header_token_b = (resolve_header_umi_token(orig_read, umi_opts, duplex_opts) or fallback_umi).encode("ascii")

                for dup_idx in range(1, repeat_count):
                    dup_qname = f"{final_qname}_boost_{dup_idx}"
                    dup_qname_b = dup_qname.encode("ascii")

                    if apply_dup_noise:
                        dup_read = r_clean.copy()
                        dup_read.query_name = dup_qname
                        dup_read.query_sequence = base_seq
                        dup_read.query_qualities = list(base_qual)
                        dup_read = inject_noise_to_record(dup_read, ref_chrom, bg_model, error_opts, mut_indices, aligned_pairs=aligned_pairs, py_rng=noise_py_rng)

                        d_seq_str = dup_read.query_sequence
                        d_qual_bytes = _quals_to_phred_bytes(dup_read.query_qualities, len(d_seq_str))
                        if dup_read.is_reverse:
                            d_seq_str = reverse_complement(d_seq_str)
                            d_qual_bytes = d_qual_bytes[::-1]

                        d_seq_b = d_seq_str.encode("ascii")
                        d_qual_b = d_qual_bytes
                    else:
                        d_seq_str = base_seq
                        d_qual_bytes = _quals_to_phred_bytes(base_qual, len(d_seq_str))
                        if r_clean.is_reverse:
                            d_seq_str = reverse_complement(d_seq_str)
                            d_qual_bytes = d_qual_bytes[::-1]

                        d_seq_b = d_seq_str.encode("ascii")
                        d_qual_b = d_qual_bytes

                    if attach_umi and placement == "inline":
                        d_seq_b = b"".join([boost_seq_prefix_b, d_seq_b])
                        d_qual_b = b"".join([boost_qual_prefix_b, d_qual_b])

                    if placement == "header":
                        header_b = b"@" + dup_qname_b + b":" + header_token_b + (b" 1:N:0:1\n" if read_num == 1 else b" 2:N:0:1\n")
                    else:
                        header_b = b"@" + dup_qname_b + (b" 1:N:0:1\n" if read_num == 1 else b" 2:N:0:1\n")

                    out_writer.write_record(b"".join([header_b, d_seq_b, b"\n+\n", d_qual_b, b"\n"]))

        orphans_flushed = 0
        pair_stats = {}
        spill_dir = os.path.dirname(os.path.abspath(out_r1_path)) or "."

        for kind, first, second in iter_mate_pairs(bam, spill_dir, threads=reader_threads, stats=pair_stats):
            if kind == "pair":
                r1_qname = strip_boost_suffix(first.query_name)
                boost_mult = read_boost_map.get(r1_qname, 1) if boost_enabled else 1

                process_and_emit_read(first, is_orphan=False, forced_boost_count=boost_mult)
                process_and_emit_read(second, is_orphan=False, forced_boost_count=boost_mult)
            else:
                orph_boost = read_boost_map.get(strip_boost_suffix(first.query_name), 1) if boost_enabled else 1
                process_and_emit_read(first, is_orphan=True, forced_boost_count=orph_boost)
                orphans_flushed += 1

        if pair_stats.get("distant_pairs_paired", 0) > 0:
            print(f"[+] Paired {pair_stats['distant_pairs_paired']:,} read pairs whose mates map to different contigs "
                  f"(or >{DISTANT_MATE_BP:,} bp apart) via a temporary name-sorted file.")
        if orphans_flushed > 0:
            print(f"[+] Total un-mated orphans routed to separate/discard stream: {orphans_flushed:,}")
    # <- writers are closed (or aborted on error) here, before the ledger is written

    if ledger_path and not no_spike:
        write_spike_ledger(ledger_path, active_variants, ledger_counts, variant_stats)

    if mutated_reads_path and not no_spike:
        with open(mutated_reads_path, "w") as f:
            for qn in sorted(mutated_qnames_set):
                f.write(f"{qn}\n")
        print(f"[+] Mutated read sidecar file written to '{mutated_reads_path}'.")

    umi_tracker.summarize()
    print("[+] Multi-SNV benchmark dataset generation complete.")
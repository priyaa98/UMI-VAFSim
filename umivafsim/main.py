# main.py
import os
import sys
import csv
import random
import numpy as np
import pysam

from .cli import parse_arguments
from .database import build_cosmic_db, fetch_cosmic_variant
from .annotator import BedRegionAnnotator
from .simulator import run_global_bam_to_fastq_multi, assert_umi_tagged
from .error_model import BackgroundNoiseModel, build_background_model, compute_provenance
from .utils import normalize_chrom

REQUIRED_MANIFEST_HEADERS = {"cosmic_id", "vaf"}


def count_manifest_rows(manifest_path):
    """Counts non-blank data rows (after the header) in a manifest; 0 if no manifest was given. Used by --no-spike."""
    if not manifest_path:
        return 0
    with open(manifest_path, "rt", encoding="utf-8") as f:
        rows = [r for r in csv.reader(f) if any(c.strip() for c in r)]
    return max(0, len(rows) - 1)


def parse_manifest(manifest_path, global_duplex_mode):
    """Parses and validates canonical multi-variant CSV manifest file."""
    print(f"[-] Parsing input variant manifest: {manifest_path}...")
    if not os.path.exists(manifest_path):
        sys.exit(f"[FATAL] Specified manifest file not found: {manifest_path}")

    variants = []
    with open(manifest_path, "rt", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = set(reader.fieldnames or [])

        missing = REQUIRED_MANIFEST_HEADERS - headers
        if missing:
            sys.exit(f"[FATAL] Manifest CSV missing required headers: {sorted(list(missing))}")

        has_spike_mode_col = "spike_mode" in headers

        if global_duplex_mode == "single" and has_spike_mode_col:
            sys.exit("[FATAL] Manifest contains 'spike_mode' column, which is invalid when --mode is 'single'.")

        # orientation_bias is duplex-only: it moves the alt between the two strands of a duplex family (/A vs /B).
        # In single-UMI data a family is one strand, so orientation can only differ BETWEEN families, and this tool
        # does not select families by orientation.
        if global_duplex_mode == "single" and "orientation_bias" in headers:
            sys.exit("[FATAL] Manifest contains an 'orientation_bias' column, which is invalid when --mode is 'single'. "
                     "orientation_bias is duplex-only: in single-UMI data each family is one strand, so there is no "
                     "within-family strand to assign. Remove the column.")

        for idx, row in enumerate(reader, start=2):
            cid = row["cosmic_id"].strip()
            if not cid:
                sys.exit(f"[FATAL] Manifest row {idx}: Missing required 'cosmic_id'.")

            try:
                vaf = float(row["vaf"].strip())
                if not (0.0 <= vaf <= 1.0):
                    raise ValueError
            except ValueError:
                sys.exit(f"[FATAL] Manifest row {idx}: Invalid 'vaf' value '{row['vaf']}'. Must be float between 0.0 and 1.0.")

            # strand_bias = fraction of the alt-carrying reads that are FORWARD-aligned (what Mutect2's SB reports).
            # Blank = no steering (every read covering the site in a mutant family carries the alt).
            sb_str = str(row.get("strand_bias", "") or "").strip()
            if not sb_str:
                strand_bias_enabled = False
                fwd_bias_ratio = 0.50
            else:
                try:
                    fwd_bias_ratio = float(sb_str)
                    if not (0.0 <= fwd_bias_ratio <= 1.0):
                        raise ValueError
                    strand_bias_enabled = True
                except ValueError:
                    sys.exit(f"[FATAL] Manifest row {idx}: Invalid 'strand_bias' ratio '{sb_str}'. Must be float in range [0.0, 1.0].")

            if has_spike_mode_col:
                raw_spike_mode = row["spike_mode"].strip().lower()
                if not raw_spike_mode:
                    sys.exit(f"[FATAL] Manifest row {idx}: 'spike_mode' column is present but cell is blank.")
                if raw_spike_mode not in ("both-strands", "single-strand"):
                    sys.exit(f"[FATAL] Manifest row {idx}: Invalid 'spike_mode' '{raw_spike_mode}'. Must be 'both-strands' or 'single-strand'.")
                spike_mode = raw_spike_mode
            else:
                spike_mode = "both-strands"

            # orientation_bias = fraction of alt-carrying families whose alt sits on strand A (the original strand
            # tagged /A), the rest on /B. Needs duplex data and single-strand spikes; blank = random 50:50 pick.
            ob_str = str(row.get("orientation_bias", "") or "").strip()
            ob_enabled, ob_val = False, 0.5
            if ob_str:
                try:
                    ob_val = float(ob_str)
                    if not (0.0 <= ob_val <= 1.0):
                        raise ValueError
                except ValueError:
                    sys.exit(f"[FATAL] Manifest row {idx}: Invalid 'orientation_bias' '{ob_str}'. Must be float in range [0.0, 1.0].")
                if global_duplex_mode != "duplex":  # unreachable in practice (column rejected above); kept as a safeguard
                    sys.exit(f"[FATAL] Manifest row {idx}: 'orientation_bias' needs --mode duplex (strand of origin comes from duplex tags).")
                if spike_mode != "single-strand":
                    sys.exit(f"[FATAL] Manifest row {idx}: 'orientation_bias' needs spike_mode 'single-strand'; with 'both-strands' the alt is on both strands of every chosen family.")
                ob_enabled = True

            variants.append({
                "cosmic_id": cid,
                "vaf": vaf,
                "strand_bias_enabled": strand_bias_enabled,
                "fwd_bias_ratio": fwd_bias_ratio,
                "spike_mode": spike_mode,
                "orientation_bias_enabled": ob_enabled,
                "orientation_bias": ob_val,
                "row_idx": idx
            })

    if not variants:
        sys.exit(f"[FATAL] Manifest file '{manifest_path}' contains no valid variant entries.")

    print(f"[+] Successfully loaded {len(variants)} target variant specification(s) from manifest.")
    return variants


def main():
    args = parse_arguments()

    random.seed(args.seed)
    np.random.seed(args.seed)

    out_dir = os.path.dirname(args.out_r1)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    if args.no_spike:
        pass  # no COSMIC lookups are needed
    elif not os.path.exists(args.db_path):
        if not args.cosmic_tar:
            sys.exit(f"[FATAL] Local database file '{args.db_path}' does not exist, and --cosmic-tar was not provided!")
        build_cosmic_db(args.cosmic_tar, args.db_path)
    else:
        print(f"[+] Found existing database store at '{args.db_path}'. Skipping build phase.")

    if args.no_spike:
        n_rows = count_manifest_rows(args.manifest)
        print("[+] --no-spike: no variants will be spiked; FASTQs are written from the unmodified BAM reads.")
        if n_rows:
            print(f"[WARNING] --no-spike: ignoring {n_rows} variant row(s) in the manifest '{args.manifest}'.")
        manifest_variants = []
    else:
        manifest_variants = parse_manifest(args.manifest, args.mode)

    with pysam.AlignmentFile(args.coord_bam, "rb") as bam:
        bam_refs = set(bam.references)

    bed_annotator = BedRegionAnnotator(args.bed_file) if (args.bed_file and not args.no_spike) else None
    resolved_variants = []
    hard_fatal_errors = []
    bed_excluded_variants = []

    for v_spec in manifest_variants:
        cid = v_spec["cosmic_id"]
        db_info = fetch_cosmic_variant(cid, v_spec["vaf"], db_path=args.db_path)

        if not db_info:
            hard_fatal_errors.append(f"Row {v_spec['row_idx']}: COSMIC ID '{cid}' not found in database.")
            continue

        raw_chrom = db_info["chrom"]
        norm_chrom = normalize_chrom(raw_chrom, bam_refs)

        if norm_chrom not in bam_refs:
            hard_fatal_errors.append(f"Row {v_spec['row_idx']}: Chromosome '{raw_chrom}' (normalized: '{norm_chrom}') missing in BAM header.")
            continue

        db_info["chrom"] = norm_chrom
        db_info.update(v_spec)

        if bed_annotator:
            is_covered, reg_tag, d_start, d_end, min_edge_d = bed_annotator.annotate_coordinate(
                norm_chrom, db_info["pos"], bam_refs=bam_refs
            )
            if not is_covered:
                bed_excluded_variants.append(db_info)
                continue

            db_info.update({
                "bed_annotated": True,
                "bed_region_tag": reg_tag,
                "dist_to_start": d_start,
                "dist_to_end": d_end,
                "min_edge_dist": min_edge_d
            })
        else:
            db_info["bed_annotated"] = False

        resolved_variants.append(db_info)

    if hard_fatal_errors:
        print("\n[FATAL] Pre-flight Manifest Validation Errors Encountered:")
        for err in hard_fatal_errors:
            print(f"  - {err}")
        sys.exit(1)

    if bed_annotator:
        if len(resolved_variants) == 0:
            sys.exit("[FATAL] All variants fall outside specified BED target regions.")
        elif bed_excluded_variants:
            print("\n[WARNING] The following variant(s) fall outside BED coverage and will be excluded:")
            for ev in bed_excluded_variants:
                print(f"  - Row {ev['row_idx']}: Variant '{ev['cosmic_id']}' ({ev['chrom']}:{ev['pos']})")
            print()

    detected_strand_source = assert_umi_tagged(
        args.coord_bam,
        resolved_variants,
        probe_window=args.umi_probe_window,
        mode=args.mode,
        umi_len=args.umi_len
    )

    bg_model = None
    if args.enable_error_model:
        control_bams = args.control_bams or [args.coord_bam]

        if args.noise_model_scope == "local":
            bed_regions = [
                (v["chrom"], max(0, v["pos"] - 1 - 500), v["pos"] - 1 + 500)
                for v in resolved_variants
            ]
        elif bed_annotator:
            bed_regions = bed_annotator.get_all_regions()
        else:
            bed_regions = [
                (v["chrom"], max(0, v["pos"] - 1 - 500), v["pos"] - 1 + 500)
                for v in resolved_variants
            ]

        if args.noise_model_cache and os.path.exists(args.noise_model_cache) and not args.rebuild_noise_model:
            bg_model = BackgroundNoiseModel.load(args.noise_model_cache)
            expected_provenance = compute_provenance(control_bams, bed_regions)
            if bg_model.provenance != expected_provenance:
                sys.exit(
                    f"[FATAL] Cached noise model at '{args.noise_model_cache}' was built from "
                    f"different control BAMs or BED regions than this run."
                )
            print(f"[+] Loaded cached background noise model from '{args.noise_model_cache}'.")
        else:
            bg_model = build_background_model(control_bams, bed_regions, seed=args.seed)
            if args.noise_model_cache:
                bg_model.save(args.noise_model_cache)
                print(f"[+] Saved background noise model to '{args.noise_model_cache}'.")

    duplex_opts = {
        "mode": args.mode,
        "strand_source": detected_strand_source,
        "boost_linkage": args.duplex_boost_linkage
    }

    orphan_opts = {
        "mode": args.orphan_mode,
        "out_orphan_r1": args.out_orphan_r1,
        "out_orphan_r2": args.out_orphan_r2
    }

    run_global_bam_to_fastq_multi(
        coord_bam=args.coord_bam,
        out_r1_path=args.out_r1,
        out_r2_path=args.out_r2,
        variants_list=resolved_variants,
        umi_opts={
            "umi_read": args.umi_read,
            "umi_len": args.umi_len,
            "spacer_len": args.spacer_len,
            "umi_placement": args.umi_placement,
            "stochastic": not args.no_stochastic,
            "qual_mode": args.qual_mode,
            "manual_qual": args.manual_qual,
            "force_noise": args.force_noise,
            "noise_threshold": args.noise_threshold,
            "background_warn_fraction": args.background_warn_fraction
        },
        boost_opts={
            "enabled": args.enable_family_boosting,
            "factor": args.boost_factor,
            "dispersion": args.boost_dispersion,
            "mode": args.boost_mode
        },
        error_opts={
            "enabled": args.enable_error_model,
            "force_global": args.force_global_error_model,
            "scale": args.error_rate_scale,
            "flat_rate": args.flat_error_rate
        },
        duplex_opts=duplex_opts,
        orphan_opts=orphan_opts,
        ledger_path=args.ledger_path,
        mutated_reads_path=args.mutated_read_names,
        bg_model=bg_model,
        bam_threads=args.threads,
        global_seed=args.seed,
        strict_strand_bias=args.strict_strand_bias,
        strict_orientation_bias=args.strict_orientation_bias,
        no_spike=args.no_spike
    )

if __name__ == "__main__":
    main()
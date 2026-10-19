# cli.py
import argparse
import os
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="In-silico multi-SNV digital spike-in engine for cfDNA and UMI-tagged sequencing datasets."
    )

    io_group = parser.add_argument_group("Required Input/Output Options")
    io_group.add_argument(
        "--coord-bam",
        required=True,
        help="Path to coordinate-sorted input BAM file containing genomic alignments."
    )
    io_group.add_argument(
        "--manifest",
        default=None,
        help="Path to CSV manifest file specifying target variants, VAFs, and injection parameters. "
             "Required unless --no-spike is given."
    )
    io_group.add_argument(
        "--out-r1",
        required=True,
        help="Path for output spiked FASTQ Read 1 file."
    )
    io_group.add_argument(
        "--out-r2",
        required=True,
        help="Path for output spiked FASTQ Read 2 file."
    )

    orphan_group = parser.add_argument_group("Orphan Output Options")
    orphan_group.add_argument(
        "--orphan-mode",
        choices=["separate", "discard"],
        default="separate",
        help="Handling of un-mated orphan reads: 'separate' writes to dedicated files, 'discard' drops them (default: separate)."
    )
    orphan_group.add_argument(
        "--out-orphan-r1",
        default=None,
        help="Path for output orphan Read 1 FASTQ (defaults to <out-r1>.orphans.fastq.gz if --orphan-mode is separate)."
    )
    orphan_group.add_argument(
        "--out-orphan-r2",
        default=None,
        help="Path for output orphan Read 2 FASTQ (defaults to <out-r2>.orphans.fastq.gz if --orphan-mode is separate)."
    )

    ledger_group = parser.add_argument_group("Spike Ledger Reporting")
    ledger_group.add_argument(
        "--ledger-path",
        default=None,
        help="Path to write CSV spike ledger report (defaults to <out-r1>.spike_ledger.csv)."
    )
    ledger_group.add_argument(
        "--mutated-read-names",
        default=None,
        help="Path to write sidecar text file containing names of mutated reads (defaults to <out-r1>.mutated_reads.txt)."
    )

    parser.add_argument(
        "--bed-file",
        default=None,
        help="Optional path to target region BED file to restrict simulation coordinates."
    )

    db_group = parser.add_argument_group("COSMIC Database Settings")
    db_group.add_argument(
        "--db-path",
        default="cosmic_local.db",
        help="Path to local SQLite database for COSMIC annotations (default: cosmic_local.db)."
    )
    db_group.add_argument(
        "--cosmic-tar",
        default=None,
        help="Path to raw COSMIC .tar archive file if initializing/building the local DB."
    )

    exec_group = parser.add_argument_group("Execution Controls")
    exec_group.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Global random seed for reproducible stochastic sampling and mutation placement (default: 42)."
    )
    exec_group.add_argument(
        "--no-stochastic",
        action="store_true",
        help="Disable stochastic binomial sampling to use deterministic read counts."
    )
    exec_group.add_argument(
        "--force-noise",
        action="store_true",
        help="Bypass script halt if estimated background noise exceeds 1%%."
    )
    exec_group.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of htslib decompression threads / compression threads (default: 4)."
    )

    exec_group.add_argument(
        "--noise-threshold",
        type=float,
        default=0.01,
        help="Read-level fraction of non-reference bases at a target site above which the locus is excluded "
             "unless --force-noise is set (default: 0.01 = 1%%). Counts every non-reference base, not only the spiked alt."
    )
    exec_group.add_argument(
        "--background-warn-fraction",
        type=float,
        default=0.2,
        help="Warn (never exclude) when unspiked reads already carrying the spiked alt exceed this fraction of the "
             "target VAF at a site (default: 0.2 = 20%%; 0 disables)."
    )

    umi_group = parser.add_argument_group("UMI & Library Architecture Parameters")
    umi_group.add_argument(
        "--umi-read",
        choices=[1, 2],
        type=int,
        default=None,
        help="Read designation containing the UMI tag - required in single mode (1 or 2)."
    )
    umi_group.add_argument(
        "--umi-len",
        type=int,
        required=True,
        help="Length of the random UMI barcode sequence in base pairs."
    )
    umi_group.add_argument(
        "--umi-placement",
        choices=["inline", "header"],
        default="inline",
        help="UMI embedding location: 'inline' prepends UMI (and optional spacer); 'header' appends to header."
    )
    umi_group.add_argument(
        "--spacer-len",
        type=int,
        default=None,
        help="Length of fixed spacer/adapter adjacent to UMI (required when --umi-placement is inline)."
    )
    umi_group.add_argument(
        "--umi-probe-window",
        type=int,
        default=500,
        help="Window (in bp) sampled by upfront UMI-tagging check (default: 500)."
    )
    umi_group.add_argument(
        "--mode",
        choices=["single", "duplex"],
        required=True,
        help="[REQUIRED] UMI library protocol type: 'single' (simplex UMI) or 'duplex' (duplex UMI)."
    )

    duplex_group = parser.add_argument_group("Duplex UMI Architecture Options")
    duplex_group.add_argument(
        "--duplex-boost-linkage",
        choices=["joint", "independent"],
        default="joint",
        help="Negative binomial family boosting behavior across strands: 'joint' or 'independent' (default: joint)."
    )

    sb_group = parser.add_argument_group("Strand-Bias Options")
    sb_group.add_argument(
        "--strict-strand-bias",
        action="store_true",
        help="Abort before writing output if any manifest 'strand_bias' request is unreachable at its site "
             "(default: warn, use the closest reachable ratio, and report it in the ledger)."
    )

    ns_group = parser.add_argument_group("No-Spike (Baseline) Options")
    ns_group.add_argument(
        "--no-spike",
        action="store_true",
        help="Write the FASTQs from the input BAM with NO variants spiked, through the same code path as a spiked run "
             "(same UMI embedding, mate pairing and orphan handling). Use it for the unspiked baseline or a negative "
             "control. The manifest is optional (missing, empty or header-only is fine); if it has variant rows they are "
             "ignored and counted in a notice, so a spiked command can be reused with one extra flag. No COSMIC database "
             "is needed, no ledger or mutated-read file is written, and --enable-error-model / --enable-family-boosting "
             "are rejected because they would alter the reads."
    )

    ob_group = parser.add_argument_group("Orientation-Bias Options")
    ob_group.add_argument(
        "--strict-orientation-bias",
        action="store_true",
        help="Abort before writing output if any manifest 'orientation_bias' request is unreachable at its site "
             "(default: warn, use the closest reachable fraction, and report it in the ledger)."
    )

    boost_group = parser.add_argument_group("Negative Binomial Family Boosting Options")
    boost_group.add_argument(
        "--enable-family-boosting",
        action="store_true",
        help="Enable Negative Binomial amplification model for overdispersed UMI family sizes."
    )
    boost_group.add_argument(
        "--boost-factor",
        type=float,
        default=2.0,
        help="Scaling multiplier for family expansion size (default: 2.0)."
    )
    boost_group.add_argument(
        "--boost-dispersion",
        type=float,
        default=0.5,
        help="Dispersion parameter size (r) for Negative Binomial family distribution (default: 0.5)."
    )
    boost_group.add_argument(
        "--boost-mode",
        choices=["all-locus", "variant-only"],
        default="all-locus",
        help="Apply family size expansion across all target locus reads or variant-bearing reads only."
    )

    error_group = parser.add_argument_group("Sequencing Error Model Options")
    error_group.add_argument(
        "--enable-error-model",
        action="store_true",
        help="Enable empirical sequencing error and substitution noise simulation."
    )
    error_group.add_argument(
        "--force-global-error-model",
        action="store_true",
        help="Apply synthetic error model noise to original reads in addition to boosted duplicates."
    )
    error_group.add_argument(
        "--error-rate-scale",
        type=float,
        default=1.0,
        help="Linear scaling factor applied to base substitution error rates (default: 1.0)."
    )
    error_group.add_argument(
        "--flat-error-rate",
        type=float,
        default=None,
        help="Apply a uniform flat background error rate overriding empirical models (e.g., 0.001)."
    )

    noise_group = parser.add_argument_group("Background Noise Model Options")
    noise_group.add_argument(
        "--control-bams",
        nargs="+",
        default=None,
        help="Control BAMs for background noise profiles - defaults to --coord-bam if omitted."
    )
    noise_group.add_argument(
        "--noise-model-cache",
        default=None,
        required=False,
        help="Path to load/save the background noise model (.json.gz)."
    )
    noise_group.add_argument(
        "--rebuild-noise-model",
        action="store_true",
        help="Force recomputation of the background noise model."
    )
    noise_group.add_argument(
        "--noise-model-scope",
        choices=["panel", "local"],
        default="panel",
        help="Scope for building noise model: 'panel' or 'local'."
    )

    qual_group = parser.add_argument_group("Phred Quality Modification Options")
    qual_group.add_argument(
        "--qual-mode",
        choices=["keep", "manual", "decay"],
        default="keep",
        help="Phred score strategy: 'keep', 'manual', or 'decay' (default: keep)."
    )
    qual_group.add_argument(
        "--manual-qual",
        type=int,
        default=30,
        help="Target Phred quality score (Q) when --qual-mode is set to 'manual' (default: 30)."
    )

    args = parser.parse_args()

    def _normalize_path(p):
        return os.path.abspath(os.path.expanduser(p)) if p else p

    args.coord_bam = _normalize_path(args.coord_bam)
    args.manifest = _normalize_path(args.manifest)
    if not args.no_spike and not args.manifest:
        parser.error("--manifest is required unless --no-spike is given.")
    args.out_r1 = _normalize_path(args.out_r1)
    args.out_r2 = _normalize_path(args.out_r2)
    args.bed_file = _normalize_path(args.bed_file)
    args.db_path = _normalize_path(args.db_path)
    args.cosmic_tar = _normalize_path(args.cosmic_tar)
    args.noise_model_cache = _normalize_path(args.noise_model_cache)
    args.out_orphan_r1 = _normalize_path(args.out_orphan_r1) if args.out_orphan_r1 else None
    args.out_orphan_r2 = _normalize_path(args.out_orphan_r2) if args.out_orphan_r2 else None
    args.ledger_path = _normalize_path(args.ledger_path) if args.ledger_path else None
    args.mutated_read_names = _normalize_path(args.mutated_read_names) if args.mutated_read_names else None

    # Default output derivative paths if unassigned
    if not args.ledger_path:
        args.ledger_path = f"{args.out_r1}.spike_ledger.csv"
    if not args.mutated_read_names:
        args.mutated_read_names = f"{args.out_r1}.mutated_reads.txt"

    if args.orphan_mode == "separate":
        if not args.out_orphan_r1:
            args.out_orphan_r1 = f"{args.out_r1}.orphans.fastq.gz"
        if not args.out_orphan_r2:
            args.out_orphan_r2 = f"{args.out_r2}.orphans.fastq.gz"

    if args.control_bams:
        args.control_bams = [_normalize_path(p) for p in args.control_bams]

    # Pre-flight early file validation
    if not os.path.exists(args.coord_bam):
        sys.exit(f"[FATAL] Specified --coord-bam file not found: {args.coord_bam}")

    if args.manifest and not os.path.exists(args.manifest):
        sys.exit(f"[FATAL] Specified --manifest file not found: {args.manifest}")

    if args.no_spike:
        if args.enable_error_model or args.enable_family_boosting:
            sys.exit("[FATAL] --no-spike writes unmodified reads, so it cannot be combined with "
                     "--enable-error-model or --enable-family-boosting.")

    if args.bed_file and not os.path.exists(args.bed_file):
        sys.exit(f"[FATAL] Specified --bed-file not found: {args.bed_file}")

    if not (0.0 < args.noise_threshold <= 1.0):
        sys.exit("[FATAL] --noise-threshold must be a fraction in (0, 1].")
    if args.background_warn_fraction < 0:
        sys.exit("[FATAL] --background-warn-fraction must be >= 0.")

    if args.enable_family_boosting and args.boost_dispersion <= 0:
        sys.exit("[FATAL] --boost-dispersion must be a positive float (> 0).")

    if args.mode == "duplex":
        if args.umi_read is not None:
            sys.exit("[FATAL] Option '--umi-read' was set, but --mode is 'duplex'.")
    else:
        if args.umi_read is None:
            sys.exit("[FATAL] Option '--umi-read' (1 or 2) is required when --mode is 'single'.")
        if args.duplex_boost_linkage != "joint":
            sys.exit("[FATAL] Option '--duplex-boost-linkage' is only valid in --mode duplex.")

    if args.umi_placement == "inline":
        if args.spacer_len is None:
            sys.exit("[FATAL] Option '--spacer-len' is required when '--umi-placement' is 'inline'.")
    elif args.umi_placement == "header":
        if args.spacer_len is not None:
            sys.exit("[FATAL] Option '--spacer-len' cannot be specified when '--umi-placement' is 'header'.")
        args.spacer_len = 0

    return args
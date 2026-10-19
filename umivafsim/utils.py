# utils.py
import sys
import math
import random
import numpy as np

DUPLEX_DELIMITER = "-"
VALID_STRAND_MAP = {
    "A": "A", "a": "A", "1": "A",
    "B": "B", "b": "B", "2": "B"
}
NUC_TO_IDX = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
BASES = ['A', 'C', 'G', 'T']

VALID_NUC_SET = set("ACGTNacgtn")

REV_COMP_TRANS = str.maketrans(
    "ACGTNacgtn",
    "TGCANtgcan"
)

PHRED_TRANS_TABLE = bytes(min(max(0, q), 60) + 33 for q in range(256))

REFSEQ_HG38_MAP = {
    "NC_000001.11": "chr1", "NC_000002.12": "chr2", "NC_000003.12": "chr3",
    "NC_000004.12": "chr4", "NC_000005.10": "chr5", "NC_000006.12": "chr6",
    "NC_000007.14": "chr7", "NC_000008.11": "chr8", "NC_000009.12": "chr9",
    "NC_000010.11": "chr10", "NC_000011.10": "chr11", "NC_000012.12": "chr12",
    "NC_000013.11": "chr13", "NC_000014.9": "chr14", "NC_000015.10": "chr15",
    "NC_000016.10": "chr16", "NC_000017.11": "chr17", "NC_000018.10": "chr18",
    "NC_000019.10": "chr19", "NC_000020.11": "chr20", "NC_000021.9": "chr21",
    "NC_000022.11": "chr22", "NC_000023.11": "chrX", "NC_000024.10": "chrY",
    "NC_012920.1": "chrM", "23": "chrX", "24": "chrY", "25": "chrM"
}


def reverse_complement(seq: str) -> str:
    """Returns canonical reverse complement of DNA sequence string."""
    return seq.translate(REV_COMP_TRANS)[::-1]


def sample_negative_binomial_multiplier(target_mean, dispersion_k, rng=None):
    if target_mean <= 1.0 or dispersion_k <= 0:
        return 1
    r = float(dispersion_k)
    p = r / (r + (target_mean - 1.0))
    sampler = rng.negative_binomial if rng is not None else np.random.negative_binomial
    sampled_boost = sampler(r, p)
    return 1 + int(sampled_boost)


def calculate_injected_qual(orig_qual_char, mode="keep", read_pos=0, read_len=150, manual_qual=30):
    mode = mode.lower()
    if mode in ("keep", "matching"):
        return orig_qual_char
    elif mode in ("manual", "fixed"):
        return chr(manual_qual + 33)
    elif mode == "decay":
        decay_factor = math.exp(-0.015 * read_pos)
        base_q = ord(orig_qual_char) - 33
        new_q = max(20, int(round(base_q * decay_factor)))
        return chr(new_q + 33)

    return orig_qual_char


def calculate_injected_qual_v2(orig_qual_char, read_pos, bg_model, rng=None):
    if bg_model is not None and hasattr(bg_model, 'cycle_qual_array') and len(bg_model.cycle_qual_array) > read_pos:
        mean_q = bg_model.cycle_qual_array[read_pos]
        if mean_q >= 0:
            rand_val = rng.normal(loc=mean_q, scale=2.5) if rng is not None else np.random.normal(loc=mean_q, scale=2.5)
            sampled_q = int(np.clip(float(rand_val), 20, 41))
            return chr(sampled_q + 33)
    return orig_qual_char


def apply_sequencing_error_model(seq, qual_str, mutated_indices=None, scale=1.0, flat_rate=0.001, py_rng=None):
    seq_list = list(seq)
    bases = ['A', 'C', 'G', 'T']
    mut_set = set(mutated_indices) if mutated_indices else set()
    rng = py_rng or random

    for i in range(len(seq_list)):
        if i in mut_set:
            continue

        q_val = ord(qual_str[i]) - 33
        error_prob = min(1.0, (10 ** (-q_val / 10.0)) * scale + flat_rate)

        if rng.random() < error_prob:
            curr_base = seq_list[i].upper()
            possible_alts = [b for b in bases if b != curr_base]
            new_base = rng.choice(possible_alts)
            seq_list[i] = new_base if seq_list[i].isupper() else new_base.lower()

    return "".join(seq_list)


def apply_sequencing_error_model_v2(
    seq, 
    qual_str, 
    chrom, 
    ref_start, 
    bg_model, 
    mutated_indices=None, 
    scale=1.0, 
    flat_rate=0.0005, 
    r_obj=None, 
    aligned_pairs=None,
    py_rng=None
):
    seq_list = list(seq)
    bases = ['A', 'C', 'G', 'T']
    mut_set = set(mutated_indices) if mutated_indices else set()
    rng = py_rng or random

    pos_to_ref = {}
    if aligned_pairs is not None:
        for q_pos, r_pos in aligned_pairs:
            if q_pos is not None and r_pos is not None:
                pos_to_ref[q_pos] = r_pos
    elif r_obj is not None and not r_obj.is_unmapped:
        try:
            for q_pos, r_pos in r_obj.get_aligned_pairs(matches_only=True):
                if q_pos is not None and r_pos is not None:
                    pos_to_ref[q_pos] = r_pos
        except ValueError:
            pass

    for i in range(len(seq_list)):
        if i in mut_set:
            continue

        genome_pos = pos_to_ref.get(i)
        curr_base = seq_list[i].upper()
        p_base = 0.0

        if genome_pos is not None and chrom is not None:
            bg_chrom_key = normalize_chrom(chrom, set(bg_model.position_probs.keys())) if (bg_model and hasattr(bg_model, 'position_probs')) else chrom
            lookup_key = (bg_chrom_key, genome_pos)
            if bg_model is not None and hasattr(bg_model, 'position_probs') and lookup_key in bg_model.position_probs:
                probs = bg_model.position_probs[lookup_key]
                if curr_base in NUC_TO_IDX:
                    match_p = probs[NUC_TO_IDX[curr_base]]
                    p_base = max(0.0, 1.0 - float(match_p))

        p_error = min(0.20, (p_base * scale) + flat_rate)

        if rng.random() < p_error:
            possible_alts = [b for b in bases if b != curr_base]
            new_base = rng.choice(possible_alts)
            seq_list[i] = new_base if seq_list[i].isupper() else new_base.lower()

    return "".join(seq_list)


def read_direction(is_reverse):
    """Aligned direction of a read: 'F' (forward) or 'R' (reverse). This is what Mutect2's SB field counts.
    (Named read_direction, not orientation: 'orientation bias' means which ORIGINAL STRAND carries the alt.)"""
    return "R" if is_reverse else "F"


read_orientation = read_direction  # backward-compatible alias


def pick_orientation(c_fwd, c_rev, alt_fwd, alt_rev, target_fwd, rng):
    """Chooses which aligned orientation ('F' or 'R') carries the alt for one family.

    c_fwd / c_rev:      reads of this family covering the site, by aligned orientation.
    alt_fwd / alt_rev:  alt reads already planned (running totals).
    target_fwd:         requested fraction of alt reads that are forward-aligned.

    A family covered by one orientation only has no choice. Otherwise the option that leaves the running
    forward fraction closest to the target wins (ties broken with rng)."""
    if c_rev == 0:
        return "F"
    if c_fwd == 0:
        return "R"
    frac_if_f = (alt_fwd + c_fwd) / (alt_fwd + alt_rev + c_fwd)
    frac_if_r = alt_fwd / (alt_fwd + alt_rev + c_rev)
    d_f, d_r = abs(frac_if_f - target_fwd), abs(frac_if_r - target_fwd)
    if abs(d_f - d_r) < 1e-12:
        return "F" if rng.random() < 0.5 else "R"
    return "F" if d_f < d_r else "R"


def strand_bias_reachable_range(fam_counts):
    """Range of forward alt-read fractions reachable for a chosen family set.

    fam_counts: iterable of (c_fwd, c_rev). Families covered by one orientation only are forced; the rest can
    go either way. Returns (min_frac, max_frac), or (None, None) when no alt reads are possible."""
    forced_f = forced_r = 0
    free_f = free_r = 0
    for c_f, c_r in fam_counts:
        if c_r == 0:
            forced_f += c_f
        elif c_f == 0:
            forced_r += c_r
        else:
            free_f += c_f
            free_r += c_r
    den_min = forced_f + forced_r + free_r
    den_max = forced_f + forced_r + free_f
    lo = (forced_f / den_min) if den_min > 0 else None
    hi = ((forced_f + free_f) / den_max) if den_max > 0 else None
    return lo, hi


def normalize_chrom(chrom, bam_references):
    if not chrom:
        return chrom
    chrom_str = str(chrom).strip()
    
    if chrom_str in bam_references:
        return chrom_str

    if chrom_str.startswith("chr") and chrom_str[3:] in bam_references:
        return chrom_str[3:]
    if not chrom_str.startswith("chr") and f"chr{chrom_str}" in bam_references:
        return f"chr{chrom_str}"

    if chrom_str in REFSEQ_HG38_MAP:
        mapped = REFSEQ_HG38_MAP[chrom_str]
        if mapped in bam_references:
            return mapped
        if mapped.startswith("chr") and mapped[3:] in bam_references:
            return mapped[3:]

    reverse_refseq_map = {v: k for k, v in REFSEQ_HG38_MAP.items()}
    if chrom_str in reverse_refseq_map and reverse_refseq_map[chrom_str] in bam_references:
        return reverse_refseq_map[chrom_str]

    return chrom_str


def parse_mi_duplex(mi_tag_val):
    if not mi_tag_val:
        return None, None
    if "/" in mi_tag_val:
        tokens = mi_tag_val.split("/", 1)
        fam_id = tokens[0]
        strand_token = tokens[1]
        strand_key = VALID_STRAND_MAP.get(strand_token)
        return fam_id, strand_key

    return mi_tag_val, None


def extract_umi_sequence(read, read_num, umi_len):
    """Extracts UMI sequence from tags (RX, MI, OX) or read header fallback."""
    if read.has_tag("RX"):
        rx_val = read.get_tag("RX")
        if rx_val:
            if DUPLEX_DELIMITER in rx_val:
                umi_str, err = split_rx_pair(rx_val, read_num, umi_len)
                if not err:
                    return umi_str
            elif len(rx_val) >= umi_len and set(rx_val[:umi_len]).issubset(VALID_NUC_SET):
                return rx_val[:umi_len].upper()

    if read.has_tag("MI"):
        mi_val = read.get_tag("MI")
        if mi_val:
            base_mi = mi_val.split("/", 1)[0]
            if len(base_mi) >= umi_len and set(base_mi[:umi_len]).issubset(VALID_NUC_SET):
                return base_mi[:umi_len].upper()

    for alt_tag in ("OX", "BX", "UM"):
        if read.has_tag(alt_tag):
            t_val = str(read.get_tag(alt_tag)).strip()
            if len(t_val) >= umi_len and set(t_val[:umi_len]).issubset(VALID_NUC_SET):
                return t_val[:umi_len].upper()

    qname = read.query_name or ""
    if ":" in qname:
        header_parts = qname.split(":")
        candidate_umi = header_parts[-1]
        if DUPLEX_DELIMITER in candidate_umi:
            parts = candidate_umi.split(DUPLEX_DELIMITER)
            sel = parts[0] if read_num == 1 else parts[-1]
            if len(sel) >= umi_len and set(sel[:umi_len]).issubset(VALID_NUC_SET):
                return sel[:umi_len].upper()
        elif len(candidate_umi) >= umi_len and set(candidate_umi[:umi_len]).issubset(VALID_NUC_SET):
            return candidate_umi[:umi_len].upper()

    return None


def split_rx_pair(rx_val, read_num, umi_len):
    if not rx_val or DUPLEX_DELIMITER not in rx_val:
        return None, "NO_DELIMITER"

    parts = rx_val.split(DUPLEX_DELIMITER)
    if len(parts) != 2:
        return None, "INVALID_FORMAT"

    u1, u2 = parts[0].strip().upper(), parts[1].strip().upper()
    if len(u1) < umi_len or len(u2) < umi_len:
        return None, "LENGTH_MISMATCH"

    if not set(u1[:umi_len]).issubset(VALID_NUC_SET) or not set(u2[:umi_len]).issubset(VALID_NUC_SET):
        return None, "INVALID_NUCLEOTIDES"

    return (u1[:umi_len] if read_num == 1 else u2[:umi_len]), None


def parse_rx_duplex(rx_val):
    if not rx_val or DUPLEX_DELIMITER not in rx_val:
        return None, None

    parts = rx_val.split(DUPLEX_DELIMITER)
    if len(parts) != 2:
        return None, None

    u1, u2 = parts[0].strip().upper(), parts[1].strip().upper()
    if not u1 or not u2:
        return None, None

    if u1 <= u2:
        fam_id = f"{u1}-{u2}"
        strand_key = "A"
    else:
        fam_id = f"{u2}-{u1}"
        strand_key = "B"

    return fam_id, strand_key


def build_fastq_header(qname, read_num, raw_umi, placement="inline"):
    if placement == "header" and raw_umi:
        return f"@{qname}:{raw_umi} {read_num}:N:0:1\n"
    return f"@{qname} {read_num}:N:0:1\n"
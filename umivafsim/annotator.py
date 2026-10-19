import os
import sys
import gzip
from collections import defaultdict
from .utils import normalize_chrom

class BedRegionAnnotator:
    """Parses and indexes target BED coordinates for fast spatial overlap querying."""
    
    def __init__(self, bed_path):
        self.bed_path = bed_path
        self.intervals = defaultdict(list)
        self._parse_and_index()

    def _parse_and_index(self):
        if not self.bed_path or not os.path.exists(self.bed_path):
            sys.exit(f"[FATAL] Specified BED file not found: {self.bed_path}")

        print(f"[-] Indexing targeted BED file: {self.bed_path}...")

        is_gz = self.bed_path.endswith(".gz")
        
        try:
            if is_gz:
                f = gzip.open(self.bed_path, "rt", encoding="utf-8")
            else:
                f = open(self.bed_path, "r", encoding="utf-8")

            with f:
                for line in f:
                    clean_line = line.strip()
                    if clean_line.startswith(("#", "track", "browser")) or not clean_line:
                        continue
                    
                    parts = clean_line.split()
                    if len(parts) < 3:
                        continue
                    
                    chrom = parts[0]
                    
                    try:
                        start, end = int(parts[1]), int(parts[2])
                    except ValueError:
                        continue
                    
                    gene_symbol = parts[3] if len(parts) >= 4 else "Target_Region"
                    annotation = f"{gene_symbol} ({chrom}:{start}-{end})"
                    self.intervals[chrom].append((start, end, annotation))
        except Exception as e:
            sys.exit(f"[FATAL] Failed to read BED file '{self.bed_path}': {e}")
                
        for chrom in self.intervals:
            self.intervals[chrom].sort(key=lambda x: (x[0], x[1]))

        total_regions = sum(len(v) for v in self.intervals.values())
        print(f"[+] Successfully loaded and indexed {total_regions:,} target intervals from BED.")

    def _resolve_chrom_key(self, chrom, bam_refs=None):
        """Flexible lookup key helper resolving against self.intervals.keys()."""
        norm_chrom = normalize_chrom(chrom, bam_refs) if bam_refs else chrom
        return normalize_chrom(norm_chrom, set(self.intervals.keys()))

    def annotate_coordinate(self, chrom, pos, bam_refs=None):
        chrom_key = self._resolve_chrom_key(chrom, bam_refs)
        
        if not chrom_key or chrom_key not in self.intervals:
            print(f"[WARNING] Chromosome '{chrom}' (resolved key: '{chrom_key}') could not be mapped to indexed BED regions.")
            return False, "Unmapped_Chromosome", -1, -1, -1

        chr_intervals = self.intervals[chrom_key]
        zero_based_pos = pos - 1
        
        for start, end, annotation in chr_intervals:
            if start <= zero_based_pos < end:
                dist_to_start = zero_based_pos - start
                dist_to_end = end - zero_based_pos
                min_dist_to_edge = min(dist_to_start, dist_to_end)
                return True, annotation, dist_to_start, dist_to_end, min_dist_to_edge

        return False, "Not_Covered", -1, -1, -1

    def get_all_regions(self):
        """Returns every indexed BED interval as flat (chrom, start, end) tuples."""
        regions = []
        for chrom, chrom_intervals in self.intervals.items():
            for start, end, _annotation in chrom_intervals:
                regions.append((chrom, start, end))
        return regions
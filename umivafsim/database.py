# database.py
import os
import sys
import tarfile
import sqlite3
import csv
import io

def build_cosmic_db(cosmic_tar_path, db_path):
    print(f"[-] Initializing SQLite database build from COSMIC archive: {cosmic_tar_path}...")
    if not os.path.exists(cosmic_tar_path):
        sys.exit(f"[FATAL] COSMIC tar archive not found: {cosmic_tar_path}")

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        CREATE TABLE IF NOT EXISTS mutations (
        cosmic_id TEXT PRIMARY KEY,
        chrom TEXT, pos INT, ref TEXT, alt TEXT,
        gene_symbol TEXT, mut_cds TEXT, mut_aa TEXT
        )
    """)

    records_inserted = 0
    with tarfile.open(cosmic_tar_path, "r:*") as tar:
        target_member = None
        for member in tar.getmembers():
            if member.name.endswith(".tsv") or member.name.endswith(".csv") or "CosmicMutantExport" in member.name:
                target_member = member
                break

        if not target_member:
            sys.exit("[FATAL] No valid TSV/CSV variant file found inside COSMIC tar archive.")

        f = tar.extractfile(target_member)
        text_stream = io.TextIOWrapper(f, encoding="utf-8", errors="replace")
        reader = csv.reader(text_stream, delimiter="\t")

        header = next(reader, None)
        if not header:
            sys.exit("[FATAL] COSMIC TSV export file is empty.")

        header_indices = {col.strip().lower(): i for i, col in enumerate(header)}

        id_idx = header_indices.get("gene_genomic_mut_id") or header_indices.get("mutation id") or header_indices.get("cosmic_id")
        legacy_idx = header_indices.get("legacy_mutation_id") or header_indices.get("legacy mutation id")
        chrom_idx = header_indices.get("chromosome") or header_indices.get("chr")
        pos_idx = header_indices.get("genomic_mutation_start") or header_indices.get("mutation start")
        ref_idx = header_indices.get("ref_seq") or header_indices.get("reference allele")
        alt_idx = header_indices.get("mut_seq") or header_indices.get("mutant allele")
        gene_idx = header_indices.get("gene_name") or header_indices.get("gene name")
        cds_idx = header_indices.get("mutation_cds") or header_indices.get("mutation cds")
        aa_idx = header_indices.get("mutation_aa") or header_indices.get("mutation aa")

        required_indices = [chrom_idx, pos_idx, ref_idx, alt_idx]
        if any(idx is None for idx in required_indices):
            sys.exit("[FATAL] Could not map required coordinate columns from COSMIC header.")

        max_required_idx = max(filter(lambda x: x is not None, [
            id_idx, legacy_idx, chrom_idx, pos_idx, ref_idx, alt_idx, gene_idx, cds_idx, aa_idx
        ]))

        batch = []
        for row in reader:
            if not row or len(row) <= max_required_idx:
                continue

            cid = row[id_idx].strip() if id_idx is not None and id_idx < len(row) else ""
            legacy = row[legacy_idx].strip() if legacy_idx is not None and legacy_idx < len(row) else ""
            chrom = row[chrom_idx].strip().replace("chr", "")
            
            try:
                pos = int(row[pos_idx].strip())
            except ValueError:
                continue

            ref = row[ref_idx].strip().upper()
            alt = row[alt_idx].strip().upper()
            gene = row[gene_idx].strip() if gene_idx is not None and gene_idx < len(row) else ""
            cds = row[cds_idx].strip() if cds_idx is not None and cds_idx < len(row) else ""
            aa = row[aa_idx].strip() if aa_idx is not None and aa_idx < len(row) else ""

            if cid or legacy:
                primary_key = cid if cid else legacy
                batch.append((primary_key, chrom, pos, ref, alt, gene, cds, aa))

            if len(batch) >= 50000:
                cursor.executemany("""
                    INSERT OR REPLACE INTO mutations
                    (cosmic_id, chrom, pos, ref, alt, gene_symbol, mut_cds, mut_aa)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, batch)
                records_inserted += len(batch)
                batch.clear()

        if batch:
            cursor.executemany("""
                    INSERT OR REPLACE INTO mutations
                    (cosmic_id, chrom, pos, ref, alt, gene_symbol, mut_cds, mut_aa)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                """, batch)
            records_inserted += len(batch)

    conn.commit()
    conn.close()
    print(f"[+] COSMIC DB build complete. Indexed {records_inserted:,} variant records in '{db_path}'.")


def fetch_cosmic_variant(cosmic_id, vaf=0.05, db_path="cosmic_local.db"):
    if not os.path.exists(db_path):
        return None

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    cursor.execute("""
        SELECT cosmic_id, chrom, pos, ref, alt, gene_symbol, mut_cds, mut_aa
        FROM mutations
        WHERE cosmic_id = ?
    """, (cosmic_id,))

    row = cursor.fetchone()
    conn.close()

    if not row:
        return None

    return {
        "cosmic_id": row[0],
        "chrom": row[1],
        "pos": row[2],
        "ref": row[3],
        "alt": row[4],
        "gene": row[5],
        "cds": row[6],
        "aa": row[7],
        "vaf": vaf
    }
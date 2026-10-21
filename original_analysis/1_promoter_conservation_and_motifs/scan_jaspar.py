#!/usr/bin/env python3
from Bio import SeqIO
from pyjaspar import jaspardb
import csv

FASTA_FILES = {
    "human": "MIF_human_hg38_TSSpm2000.fa",
    "chimp": "MIF_chimp_panTro6_TSSpm2000.fa",
    "mouse": "MIF_mouse_mm10_TSSpm2000.fa",
}
THRESHOLD_PCT = 0.90
OUTPUT_CSV = "jaspar_scan_results.csv"

jdb = jaspardb(release="JASPAR2024")
motif_list = jdb.fetch_motifs(collection="CORE", tax_group=["vertebrates"])
print(f"Loaded {len(motif_list)} JASPAR CORE vertebrate motifs")

results = []
for species, path in FASTA_FILES.items():
    record = next(SeqIO.parse(path, "fasta"))
    seq = record.seq.upper()
    print(f"Scanning {species} ({len(seq)} bp)...")

    for m in motif_list:
        pwm = m.counts.normalize(pseudocounts=0.8)
        pssm = pwm.log_odds()
        abs_threshold = pssm.min + THRESHOLD_PCT * (pssm.max - pssm.min)
        for position, score in pssm.search(seq, threshold=abs_threshold):
            strand = "+" if position >= 0 else "-"
            pos = abs(position)
            results.append([species, m.name, m.matrix_id, pos, strand, round(score, 2)])

with open(OUTPUT_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["species", "TF_name", "matrix_id", "position", "strand", "score"])
    writer.writerows(results)

print(f"Done: {len(results)} hits written to {OUTPUT_CSV}")

from Bio import SeqIO
from Bio.Seq import Seq

mirnas = list(SeqIO.parse("candidate_miRNA_sequences_verified.fa", "fasta"))
lncrnas = list(SeqIO.parse("candidate_lncRNA_sequences.fa", "fasta"))

print(f"Loaded {len(mirnas)} miRNAs and {len(lncrnas)} lncRNAs\n")

def seed_match(mirna_seq, target_seq):
    # Seed = positions 2-8 of the mature miRNA (standard canonical seed region)
    seed_rna = str(mirna_seq[1:8]).upper()
    seed_dna = seed_rna.replace("U", "T")
    seed_rc = str(Seq(seed_dna).reverse_complement())
    target_str = str(target_seq).upper()

    positions = []
    start = 0
    while True:
        idx = target_str.find(seed_rc, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions

results = []
for mir in mirnas:
    for lnc in lncrnas:
        hits = seed_match(mir.seq, lnc.seq)
        if hits:
            results.append((mir.id, lnc.id, len(hits), hits))
            print(f"{mir.id} <-> {lnc.id}: {len(hits)} seed match(es) at position(s) {hits}")

if not results:
    print("No 7-mer seed matches found for any pair.")

with open("seed_match_results.csv", "w") as f:
    f.write("miRNA,lncRNA,n_matches,positions\n")
    for mir_id, lnc_id, n, pos in results:
        f.write(f"{mir_id},{lnc_id},{n},\"{pos}\"\n")

print(f"\nSaved {len(results)} matching pairs to seed_match_results.csv")

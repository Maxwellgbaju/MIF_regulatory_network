import json

with open("figure_mirna_lncrna_results.json") as f:
    raw = json.load(f)

best_hits = {}
all_hits_by_mirna = {}

for mirna, lines in raw.items():
    if not lines:
        best_hits[mirna] = None
        continue
    parsed = []
    for line in lines:
        if line.startswith("#") or line.startswith("miRNAid") or "haven't been set" in line:
            continue
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        try:
            gene_name = parts[3]
            gene_type = parts[4]
            clip_num = int(parts[9])
            phylop = float(parts[16]) if parts[16] not in ("", "NA") else None
        except (IndexError, ValueError):
            continue
        if gene_type != "lncRNA":
            continue
        parsed.append({"lncRNA": gene_name, "clipExpNum": clip_num, "phyloP": phylop})
    all_hits_by_mirna[mirna] = parsed
    if parsed:
        parsed.sort(key=lambda x: x["clipExpNum"], reverse=True)
        best_hits[mirna] = parsed[0]
    else:
        best_hits[mirna] = None

print("=== Best-supported lncRNA target per miRNA (by CLIP-seq experiment count) ===")
for mirna, hit in best_hits.items():
    if hit:
        print(f"  {mirna:20s} -> {hit['lncRNA']:15s} (clipExpNum={hit['clipExpNum']}, phyloP={hit['phyloP']})")
    else:
        print(f"  {mirna:20s} -> NO ENCORI LNCRNA COVERAGE (query failed or returned no lncRNA-type hits)")

# Check for convergent lncRNA targets shared across multiple miRNAs
print("\n=== Checking for shared lncRNA targets across miRNAs (top-3 per miRNA) ===")
lncrna_to_mirnas = {}
for mirna, hits in all_hits_by_mirna.items():
    for h in sorted(hits, key=lambda x: x["clipExpNum"], reverse=True)[:3]:
        lncrna_to_mirnas.setdefault(h["lncRNA"], []).append(mirna)

shared = {k: v for k, v in lncrna_to_mirnas.items() if len(set(v)) > 1}
for lncrna, mirnas in shared.items():
    print(f"  {lncrna}: targeted by {sorted(set(mirnas))}")

with open("best_lncrna_hits.json", "w") as f:
    json.dump(best_hits, f, indent=2)
print("\nSaved best_lncrna_hits.json")

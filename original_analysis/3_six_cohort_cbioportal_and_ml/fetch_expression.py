import requests
import json
import csv

studies = [
    "kirc_tcga_pan_can_atlas_2018",
    "lihc_tcga_pan_can_atlas_2018",
    "stad_tcga_pan_can_atlas_2018",
    "brca_tcga_pan_can_atlas_2018",
    "lusc_tcga_pan_can_atlas_2018",
    "skcm_tcga_pan_can_atlas_2018",
]

# Entrez Gene IDs confirmed: MIF=4282. SP1, HIF1A, RELA below are commonly-cited
# IDs but NOT yet verified in this conversation -- verify before trusting results.
genes = {
    "MIF": 4282,
    "SP1": 6667,
    "HIF1A": 3091,
    "RELA": 5970,
}

all_rows = []
for study in studies:
    profile_id = f"{study}_rna_seq_v2_mrna_median_Zscores"
    sample_list_id = f"{study}_all"
    url = f"https://www.cbioportal.org/api/molecular-profiles/{profile_id}/molecular-data/fetch?projection=SUMMARY"
    body = {
        "entrezGeneIds": list(genes.values()),
        "sampleListId": sample_list_id,
    }
    print(f"Fetching {study} ...")
    r = requests.post(url, json=body, timeout=30)
    if r.status_code != 200:
        print(f"  FAILED: HTTP {r.status_code} - {r.text[:200]}")
        continue
    data = r.json()
    print(f"  -> {len(data)} rows")
    all_rows.extend(data)

# Reverse-map entrez ID back to gene symbol for readability
id_to_symbol = {v: k for k, v in genes.items()}
for row in all_rows:
    row["geneSymbol"] = id_to_symbol.get(row["entrezGeneId"], "UNKNOWN")

with open("MIF_axis_expression_6cancers.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["studyId", "geneSymbol", "entrezGeneId", "sampleId", "patientId", "value"])
    writer.writeheader()
    for row in all_rows:
        writer.writerow({k: row.get(k) for k in writer.fieldnames})

print(f"\nSaved {len(all_rows)} total rows to MIF_axis_expression_6cancers.csv")

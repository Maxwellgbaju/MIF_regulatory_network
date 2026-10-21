import requests, csv

studies = ["kirc_tcga_pan_can_atlas_2018","lihc_tcga_pan_can_atlas_2018",
           "stad_tcga_pan_can_atlas_2018","brca_tcga_pan_can_atlas_2018",
           "lusc_tcga_pan_can_atlas_2018","skcm_tcga_pan_can_atlas_2018"]

genes = {"MIF": 4282, "SP1": 6667, "HIF1A": 3091, "RELA": 5970}

rows = []
for study in studies:
    profile_id = f"{study}_mutations"
    sample_list_id = f"{study}_all"
    url = f"https://www.cbioportal.org/api/molecular-profiles/{profile_id}/mutations/fetch?projection=SUMMARY"
    body = {"entrezGeneIds": list(genes.values()), "sampleListId": sample_list_id}
    r = requests.post(url, json=body, timeout=30)
    if r.status_code != 200:
        print(f"{study}: FAILED HTTP {r.status_code}")
        continue
    data = r.json()
    print(f"{study}: {len(data)} mutation records")
    id_to_symbol = {v: k for k, v in genes.items()}
    for row in data:
        rows.append({"studyId": row["studyId"], "geneSymbol": id_to_symbol[row["entrezGeneId"]],
                     "sampleId": row["sampleId"], "proteinChange": row.get("proteinChange", ""),
                     "mutationType": row.get("mutationType", "")})

with open("mutations_6cancers.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["studyId","geneSymbol","sampleId","proteinChange","mutationType"])
    w.writeheader()
    w.writerows(rows)
print(f"\nSaved {len(rows)} rows")

# Also need total sample count per study to compute mutation frequency (%)
print("\nFetching total sample counts per study for frequency calculation...")
sample_counts = {}
for study in studies:
    url = f"https://www.cbioportal.org/api/studies/{study}/sample-lists/{study}_all"
    r = requests.get(url, timeout=30)
    if r.status_code == 200:
        sample_counts[study] = len(r.json()["sampleIds"])
        print(f"{study}: {sample_counts[study]} total samples")

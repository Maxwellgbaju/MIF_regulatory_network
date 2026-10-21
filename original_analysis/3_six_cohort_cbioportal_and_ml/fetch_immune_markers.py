import requests, csv

studies = ["kirc_tcga_pan_can_atlas_2018","lihc_tcga_pan_can_atlas_2018",
           "stad_tcga_pan_can_atlas_2018","brca_tcga_pan_can_atlas_2018",
           "lusc_tcga_pan_can_atlas_2018","skcm_tcga_pan_can_atlas_2018"]

# Confirmed Entrez IDs: PDCD1 (exhausted CD8+ T cell marker) = 5133,
# FOXP3 (Treg marker) = 50943, CD163 (M2 macrophage marker) = 9332
markers = {"PDCD1": 5133, "FOXP3": 50943, "CD163": 9332}

rows = []
for study in studies:
    profile_id = f"{study}_rna_seq_v2_mrna_median_Zscores"
    sample_list_id = f"{study}_all"
    url = f"https://www.cbioportal.org/api/molecular-profiles/{profile_id}/molecular-data/fetch?projection=SUMMARY"
    body = {"entrezGeneIds": list(markers.values()), "sampleListId": sample_list_id}
    print(f"Fetching {study} ...")
    r = requests.post(url, json=body, timeout=30)
    if r.status_code != 200:
        print(f"  FAILED: HTTP {r.status_code}")
        continue
    data = r.json()
    print(f"  -> {len(data)} rows")
    id_to_symbol = {v: k for k, v in markers.items()}
    for row in data:
        rows.append({"studyId": row["studyId"], "geneSymbol": id_to_symbol[row["entrezGeneId"]],
                     "sampleId": row["sampleId"], "patientId": row["patientId"], "value": row["value"]})

with open("immune_markers_6cancers.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=["studyId","geneSymbol","sampleId","patientId","value"])
    w.writeheader()
    w.writerows(rows)
print(f"\nSaved {len(rows)} rows")

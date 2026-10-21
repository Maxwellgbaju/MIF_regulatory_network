import requests
import pandas as pd
import numpy as np

genes = {
    "MIF": 4282, "RELA": 5970, "EGR1": 1958, "SP1": 6667, "NR1H2": 7376,
    "RXRA": 6256, "IRF1": 3659, "KLF4": 9314, "IRF2": 3660, "NFIL3": 4783,
    "ESR1": 2099, "TP53": 7157, "TFAP2A": 7020, "NORAD": 647979
}

study = "kirc_tcga_pan_can_atlas_2018"
profile_id = f"{study}_rna_seq_v2_mrna_median_Zscores"
sample_list_id = f"{study}_all"

print(f"Fetching {len(genes)} genes for KIRC (cBioPortal)...")
url = f"https://www.cbioportal.org/api/molecular-profiles/{profile_id}/molecular-data/fetch?projection=SUMMARY"
body = {"entrezGeneIds": list(genes.values()), "sampleListId": sample_list_id}
r = requests.post(url, json=body, timeout=30)
data = r.json()
print(f"Retrieved {len(data)} rows")

id_to_symbol = {v: k for k, v in genes.items()}
rows = []
for row in data:
    rows.append({"sampleId": row["sampleId"], "patientId": row["patientId"],
                 "gene": id_to_symbol[row["entrezGeneId"]], "value": row["value"]})
expr_df = pd.DataFrame(rows)
expr_wide = expr_df.pivot_table(index=["sampleId","patientId"], columns="gene", values="value").reset_index()
print(f"Wide expression matrix shape: {expr_wide.shape}")
expr_wide.to_csv("capstone_cbioportal_kirc_expr.csv", index=False)

# Reuse existing survival data
surv_raw = pd.read_csv("survival_6cancers_raw.csv")
surv = surv_raw.pivot_table(index=["studyId","patientId"], columns="attribute", values="value", aggfunc="first").reset_index()
surv["event"] = surv["OS_STATUS"].astype(str).str.startswith("1").astype(int)
surv["time"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")
surv = surv[surv["studyId"] == study].dropna(subset=["time"])

merged = expr_wide.merge(surv, on="patientId", how="inner")
gene_cols = list(genes.keys())
merged_clean = merged.dropna(subset=gene_cols + ["time", "event"])
print(f"\nComplete cBioPortal KIRC dataset (14 genes + survival, no missing): n={len(merged_clean)}, events={int(merged_clean['event'].sum())}")
merged_clean.to_csv("capstone_cbioportal_kirc_full.csv", index=False)
print("Saved capstone_cbioportal_kirc_full.csv")

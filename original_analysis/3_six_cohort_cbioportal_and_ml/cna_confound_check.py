import pandas as pd
from scipy.stats import spearmanr

cna = pd.read_csv("cna_6cancers.csv")
expr = pd.read_csv("MIF_axis_expression_6cancers.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}

cna_wide = cna.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
merged = cna_wide.merge(expr_wide, on=["studyId","sampleId"], suffixes=("_cna","_expr"))
print(f"Merged shape: {merged.shape}")

genes = ["MIF", "SP1", "HIF1A", "RELA"]

# --- Check A: validity check -- does each gene's own CNA status track its own expression? ---
print("\n=== A. Own-gene CNA status vs. own-gene expression (validity check) ===")
for study in merged["studyId"].unique():
    for gene in genes:
        sub = merged[[f"{gene}_cna", f"{gene}_expr"]].loc[merged["studyId"] == study].dropna()
        r, p = spearmanr(sub[f"{gene}_cna"], sub[f"{gene}_expr"])
        print(f"{study_to_code[study]:6s} {gene:6s}: rho={r:+.3f}  p={p:.2e}  (n={len(sub)})")

# --- Check B: do MIF and SP1 CNA statuses co-occur? (tests shared-CNA-event explanation) ---
print("\n=== B. MIF CNA status vs. SP1 CNA status (co-occurrence check) ===")
for study in merged["studyId"].unique():
    sub = merged[["MIF_cna", "SP1_cna"]].loc[merged["studyId"] == study].dropna()
    r, p = spearmanr(sub["MIF_cna"], sub["SP1_cna"])
    print(f"{study_to_code[study]:6s}: rho={r:+.3f}  p={p:.2e}  (n={len(sub)})")

# --- Check C: does the MIF-SP1 EXPRESSION correlation change when restricting to CNA-diploid samples only? ---
print("\n=== C. MIF-SP1 expression correlation: all samples vs. CNA-diploid-only samples ===")
for study in merged["studyId"].unique():
    all_sub = merged[["MIF_expr", "SP1_expr"]].loc[merged["studyId"] == study].dropna()
    r_all, p_all = spearmanr(all_sub["MIF_expr"], all_sub["SP1_expr"])

    diploid = merged.loc[(merged["studyId"] == study) & (merged["MIF_cna"] == 0) & (merged["SP1_cna"] == 0)]
    diploid = diploid[["MIF_expr", "SP1_expr"]].dropna()
    if len(diploid) >= 10:
        r_dip, p_dip = spearmanr(diploid["MIF_expr"], diploid["SP1_expr"])
        print(f"{study_to_code[study]:6s}: all rho={r_all:+.3f} (n={len(all_sub)})  |  diploid-only rho={r_dip:+.3f} (n={len(diploid)})")
    else:
        print(f"{study_to_code[study]:6s}: all rho={r_all:+.3f} (n={len(all_sub)})  |  diploid-only: insufficient n ({len(diploid)})")

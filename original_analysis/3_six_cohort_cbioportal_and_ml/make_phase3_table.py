import pandas as pd

corr = pd.read_csv("MIF_axis_correlations_6cancers.csv")
partial = pd.read_csv("SP1_MIF_partial_correlation.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC",
    "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD",
    "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC",
    "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
corr["Cancer"] = corr["study"].map(study_to_code)
partial["Cancer"] = partial["study"].map(study_to_code)

# --- Table 1: full 18-test correlation table (all genes x all cancers) ---
t1 = corr[["Cancer", "gene", "n", "rho", "pval", "fdr"]].copy()
t1.columns = ["Cancer", "Gene", "n", "Spearman rho", "p-value", "FDR (BH)"]
t1["Spearman rho"] = t1["Spearman rho"].round(3)
t1["p-value"] = t1["p-value"].apply(lambda x: f"{x:.2e}")
t1["FDR (BH)"] = t1["FDR (BH)"].apply(lambda x: f"{x:.2e}")
t1 = t1.sort_values(["Gene", "Cancer"]).reset_index(drop=True)

# --- Table 2: SP1-MIF raw vs. purity-adjusted partial correlation ---
t2 = partial[["Cancer", "n", "raw_rho", "raw_p", "partial_rho", "partial_p"]].copy()
t2.columns = ["Cancer", "n", "Raw rho", "Raw p", "Partial rho\n(purity-adjusted)", "Partial p"]
t2["Raw rho"] = t2["Raw rho"].round(3)
t2["Partial rho\n(purity-adjusted)"] = t2["Partial rho\n(purity-adjusted)"].round(3)
t2["Raw p"] = t2["Raw p"].apply(lambda x: f"{x:.2e}")
t2["Partial p"] = t2["Partial p"].apply(lambda x: f"{x:.2e}")
t2 = t2.sort_values("Cancer").reset_index(drop=True)

# --- Save both as CSV and as a combined Excel workbook (one tab per table) ---
t1.to_csv("Table1_MIF_axis_correlations_formatted.csv", index=False)
t2.to_csv("Table2_SP1_purity_adjustment_formatted.csv", index=False)

with pd.ExcelWriter("Phase3_Tables.xlsx", engine="openpyxl") as writer:
    t1.to_excel(writer, sheet_name="Table1_Correlations", index=False)
    t2.to_excel(writer, sheet_name="Table2_PuritySensitivity", index=False)

print("=== Table 1: MIF-axis gene correlations across 6 cancers ===")
print(t1.to_string(index=False))
print("\n=== Table 2: SP1-MIF correlation, raw vs. purity-adjusted ===")
print(t2.to_string(index=False))
print("\nSaved: Table1_MIF_axis_correlations_formatted.csv, Table2_SP1_purity_adjustment_formatted.csv, Phase3_Tables.xlsx")

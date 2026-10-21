import pandas as pd
import numpy as np

corr = pd.read_csv("MIF_axis_correlations_6cancers.csv")
partial = pd.read_csv("SP1_MIF_partial_correlation.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
corr["Cancer"] = corr["study"].map(study_to_code)
partial["Cancer"] = partial["study"].map(study_to_code)

def fisher_ci(rho, n, alpha=0.05):
    z = np.arctanh(rho)
    se = 1 / np.sqrt(n - 3)
    zcrit = 1.959963984540054
    lo, hi = np.tanh(z - zcrit * se), np.tanh(z + zcrit * se)
    return lo, hi

# --- Table 1 ---
corr["ci_lo"], corr["ci_hi"] = zip(*corr.apply(lambda r: fisher_ci(r["rho"], r["n"]), axis=1))
t1 = corr[["Cancer", "gene", "n", "rho", "ci_lo", "ci_hi", "pval", "fdr"]].copy()
t1.columns = ["Cancer", "Gene", "n", "Spearman rho", "95% CI lower", "95% CI upper", "p-value", "FDR (BH)"]
for c in ["Spearman rho", "95% CI lower", "95% CI upper"]:
    t1[c] = t1[c].round(3)
t1["p-value"] = t1["p-value"].apply(lambda x: f"{x:.2e}")
t1["FDR (BH)"] = t1["FDR (BH)"].apply(lambda x: f"{x:.2e}")
t1 = t1.sort_values(["Gene", "Cancer"]).reset_index(drop=True)

# --- Table 2 ---
partial["raw_ci_lo"], partial["raw_ci_hi"] = zip(*partial.apply(lambda r: fisher_ci(r["raw_rho"], r["n"]), axis=1))
partial["partial_ci_lo"], partial["partial_ci_hi"] = zip(*partial.apply(lambda r: fisher_ci(r["partial_rho"], r["n"]), axis=1))
t2 = partial[["Cancer", "n", "raw_rho", "raw_ci_lo", "raw_ci_hi", "raw_p",
              "partial_rho", "partial_ci_lo", "partial_ci_hi", "partial_p"]].copy()
t2.columns = ["Cancer", "n", "Raw rho", "Raw 95% CI lower", "Raw 95% CI upper", "Raw p",
              "Partial rho (purity-adjusted)", "Partial 95% CI lower", "Partial 95% CI upper", "Partial p"]
for c in ["Raw rho","Raw 95% CI lower","Raw 95% CI upper","Partial rho (purity-adjusted)","Partial 95% CI lower","Partial 95% CI upper"]:
    t2[c] = t2[c].round(3)
t2["Raw p"] = t2["Raw p"].apply(lambda x: f"{x:.2e}")
t2["Partial p"] = t2["Partial p"].apply(lambda x: f"{x:.2e}")
t2 = t2.sort_values("Cancer").reset_index(drop=True)

t1.to_csv("Table1_MIF_axis_correlations_formatted.csv", index=False)
t2.to_csv("Table2_SP1_purity_adjustment_formatted.csv", index=False)

with pd.ExcelWriter("Phase3_Tables.xlsx", engine="openpyxl") as writer:
    t1.to_excel(writer, sheet_name="Table1_Correlations", index=False)
    t2.to_excel(writer, sheet_name="Table2_PuritySensitivity", index=False)

print("=== Table 1 ===")
print(t1.to_string(index=False))
print("\n=== Table 2 ===")
print(t2.to_string(index=False))
print("\nSaved corrected files.")

import pandas as pd
import numpy as np

res = pd.read_csv("MIF_axis_immune_correlations_partial.csv")
boot = pd.read_csv("immune_bootstrap_diff_all.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
res["Cancer"] = res["study"].map(study_to_code)
boot["Cancer"] = boot["study"].map(study_to_code)

merged = res.merge(boot[["study","axis_gene","immune_marker","diff","diff_ci_lo","diff_ci_hi","bootstrap_p"]],
                    on=["study","axis_gene","immune_marker"], how="left")

def fmt_p(x):
    return f"{x:.2e}" if x < 0.001 else f"{x:.3f}"

# --- Table A: PDCD1 headline results (main text) ---
main = merged[merged["immune_marker"] == "PDCD1"].copy()
main = main[["Cancer","axis_gene","n","raw_rho","partial_rho_ctrl_PTPRC","partial_fdr","diff","diff_ci_lo","diff_ci_hi","bootstrap_p"]]
main.columns = ["Cancer","Gene","n","Raw rho","Partial rho (ctrl. PTPRC)","Partial FDR","Raw-Partial diff","Diff 95% CI lower","Diff 95% CI upper","Bootstrap p (shift)"]
for c in ["Raw rho","Partial rho (ctrl. PTPRC)","Raw-Partial diff","Diff 95% CI lower","Diff 95% CI upper"]:
    main[c] = main[c].round(3)
main["Partial FDR"] = main["Partial FDR"].apply(fmt_p)
main["Bootstrap p (shift)"] = main["Bootstrap p (shift)"].apply(fmt_p)
cancer_order = ["KIRC","LIHC","STAD","BRCA","LUSC","SKCM"]
gene_order = ["MIF","SP1","HIF1A","RELA"]
main["Cancer"] = pd.Categorical(main["Cancer"], categories=cancer_order, ordered=True)
main["Gene"] = pd.Categorical(main["Gene"], categories=gene_order, ordered=True)
main = main.sort_values(["Gene","Cancer"]).reset_index(drop=True)
main.to_csv("Table_PDCD1_main.csv", index=False)

# --- Table B: full 72-test supplementary table ---
supp = merged[["Cancer","axis_gene","immune_marker","n","raw_rho","partial_rho_ctrl_PTPRC","partial_fdr","diff","bootstrap_p"]].copy()
supp.columns = ["Cancer","Gene","Immune marker","n","Raw rho","Partial rho (ctrl. PTPRC)","Partial FDR","Raw-Partial diff","Bootstrap p (shift)"]
for c in ["Raw rho","Partial rho (ctrl. PTPRC)","Raw-Partial diff"]:
    supp[c] = supp[c].round(3)
supp["Partial FDR"] = supp["Partial FDR"].apply(fmt_p)
supp["Bootstrap p (shift)"] = supp["Bootstrap p (shift)"].apply(fmt_p)
supp["Cancer"] = pd.Categorical(supp["Cancer"], categories=cancer_order, ordered=True)
supp["Gene"] = pd.Categorical(supp["Gene"], categories=gene_order, ordered=True)
supp["Immune marker"] = pd.Categorical(supp["Immune marker"], categories=["PDCD1","FOXP3","CD163"], ordered=True)
supp = supp.sort_values(["Immune marker","Gene","Cancer"]).reset_index(drop=True)
supp.to_csv("Table_all72_supplementary.csv", index=False)

with pd.ExcelWriter("Phase3_Step3.3_Tables.xlsx", engine="openpyxl") as writer:
    main.to_excel(writer, sheet_name="Main_PDCD1", index=False)
    supp.to_excel(writer, sheet_name="Supplementary_All72", index=False)

print("=== Main table (PDCD1, 24 rows) ===")
print(main.to_string(index=False))
print(f"\nSaved: Table_PDCD1_main.csv, Table_all72_supplementary.csv, Phase3_Step3.3_Tables.xlsx")

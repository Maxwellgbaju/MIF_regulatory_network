import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

cna_freq = pd.read_csv("CNA_frequency_summary.csv")
mut_freq = pd.read_csv("Mutation_frequency_summary.csv")
cna = pd.read_csv("cna_6cancers.csv")
expr = pd.read_csv("MIF_axis_expression_6cancers.csv")

cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["MIF", "SP1", "HIF1A", "RELA"]

def fisher_ci(rho, n):
    z = np.arctanh(rho)
    se = 1 / np.sqrt(n - 3)
    zc = 1.959963984540054
    return np.tanh(z - zc*se), np.tanh(z + zc*se)

fig, axes = plt.subplots(1, 2, figsize=(15, 6), gridspec_kw={"width_ratios": [1.3, 1]})

# ---- Panel A: CNA gain/loss frequency heatmap (net alteration = gain% - loss%) ----
ax0 = axes[0]
cna_freq["net"] = cna_freq["pct_any_gain"] - cna_freq["pct_any_loss"]
mat = cna_freq.pivot(index="Cancer", columns="Gene", values="net").reindex(index=cancer_order, columns=gene_order)
gain_mat = cna_freq.pivot(index="Cancer", columns="Gene", values="pct_any_gain").reindex(index=cancer_order, columns=gene_order)
loss_mat = cna_freq.pivot(index="Cancer", columns="Gene", values="pct_any_loss").reindex(index=cancer_order, columns=gene_order)
mut_mat = mut_freq.pivot(index="Cancer", columns="Gene", values="pct_mutated").reindex(index=cancer_order, columns=gene_order)

vmax = np.nanmax(np.abs(mat.values))
im = ax0.imshow(mat.values, cmap="PuOr_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax0.set_xticks(range(len(gene_order))); ax0.set_xticklabels(gene_order, fontsize=11)
ax0.set_yticks(range(len(cancer_order))); ax0.set_yticklabels(cancer_order, fontsize=11)
for i in range(len(cancer_order)):
    for j in range(len(gene_order)):
        g, l, m = gain_mat.values[i,j], loss_mat.values[i,j], mut_mat.values[i,j]
        color = "white" if abs(mat.values[i,j]) > vmax*0.5 else "black"
        ax0.text(j, i, f"+{g:.0f}%/-{l:.0f}%\nmut {m:.1f}%", ha="center", va="center", fontsize=8, color=color)
cbar = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)
cbar.set_label("Net CNA (%gain − %loss)", fontsize=9)
ax0.set_title("A. Copy-number and mutation frequency", fontsize=12, loc="left")

# ---- Panel B: SP1-MIF correlation, all samples vs. CNA-diploid-only, with 95% CI ----
cna_wide = cna.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
merged = cna_wide.merge(expr_wide, on=["studyId","sampleId"], suffixes=("_cna","_expr"))
study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}

results = []
for study, code in study_to_code.items():
    all_sub = merged[["MIF_expr","SP1_expr"]].loc[merged["studyId"]==study].dropna()
    r_all, _ = spearmanr(all_sub["MIF_expr"], all_sub["SP1_expr"])
    dip = merged.loc[(merged["studyId"]==study) & (merged["MIF_cna"]==0) & (merged["SP1_cna"]==0), ["MIF_expr","SP1_expr"]].dropna()
    r_dip, _ = spearmanr(dip["MIF_expr"], dip["SP1_expr"])
    results.append({"cancer": code, "n_all": len(all_sub), "rho_all": r_all, "n_dip": len(dip), "rho_dip": r_dip})

res_df = pd.DataFrame(results).set_index("cancer").reindex(cancer_order)

ax1 = axes[1]
x = np.arange(len(cancer_order)); width = 0.35
all_ci = [fisher_ci(res_df.loc[c,"rho_all"], res_df.loc[c,"n_all"]) for c in cancer_order]
dip_ci = [fisher_ci(res_df.loc[c,"rho_dip"], res_df.loc[c,"n_dip"]) for c in cancer_order]
all_vals = res_df["rho_all"].values
dip_vals = res_df["rho_dip"].values
all_err = [[all_vals[i]-all_ci[i][0] for i in range(6)],[all_ci[i][1]-all_vals[i] for i in range(6)]]
dip_err = [[dip_vals[i]-dip_ci[i][0] for i in range(6)],[dip_ci[i][1]-dip_vals[i] for i in range(6)]]

ax1.bar(x-width/2, all_vals, width, yerr=all_err, capsize=4, label="All samples", color="#C00000", error_kw={"elinewidth":1.2})
ax1.bar(x+width/2, dip_vals, width, yerr=dip_err, capsize=4, label="CNA-diploid only\n(MIF & SP1)", color="#2E7D32", error_kw={"elinewidth":1.2})
ax1.axhline(0, color="black", linewidth=0.8)
ax1.set_xticks(x); ax1.set_xticklabels(cancer_order, fontsize=10)
ax1.set_ylabel("Spearman ρ (SP1 vs MIF), 95% CI", fontsize=10)
ax1.set_title("B. SP1-MIF correlation persists without CNA", fontsize=12, loc="left")
ax1.legend(fontsize=9, loc="lower right")

plt.tight_layout()
plt.savefig("Phase3_4_CNA_mutation_figure.png", dpi=300, bbox_inches="tight")
plt.savefig("Phase3_4_CNA_mutation_figure.eps", format="eps", bbox_inches="tight")
print("Saved Phase3_4_CNA_mutation_figure.png/.eps")

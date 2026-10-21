import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- Load your actual results (no hardcoded numbers) ----
corr = pd.read_csv("MIF_axis_correlations_6cancers.csv")
partial = pd.read_csv("SP1_MIF_partial_correlation.csv")
boot = pd.read_csv("SP1_MIF_bootstrap_diff_test.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
corr["cancer"] = corr["study"].map(study_to_code)
partial["cancer"] = partial["study"].map(study_to_code)
boot["cancer"] = boot["study"].map(study_to_code)

cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["SP1", "HIF1A", "RELA"]

def fisher_ci(rho, n):
    z = np.arctanh(rho)
    se = 1 / np.sqrt(n - 3)
    zc = 1.959963984540054
    return np.tanh(z - zc*se), np.tanh(z + zc*se)

def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "n.s."

# ---- Panel A data ----
rho_matrix = corr.pivot(index="cancer", columns="gene", values="rho").reindex(index=cancer_order, columns=gene_order)
fdr_matrix = corr.pivot(index="cancer", columns="gene", values="fdr").reindex(index=cancer_order, columns=gene_order)

# ---- Panel B data ----
partial_ordered = partial.set_index("cancer").reindex(cancer_order)
boot_ordered = boot.set_index("cancer").reindex(cancer_order)

fig, axes = plt.subplots(1, 2, figsize=(13.5, 6), gridspec_kw={"width_ratios": [1.1, 1]})

# ---- Panel A: heatmap ----
ax0 = axes[0]
vmax = np.nanmax(np.abs(rho_matrix.values))
im = ax0.imshow(rho_matrix.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax0.set_xticks(range(len(gene_order))); ax0.set_xticklabels(gene_order, fontsize=11)
ax0.set_yticks(range(len(cancer_order))); ax0.set_yticklabels(cancer_order, fontsize=11)
for i in range(len(cancer_order)):
    for j in range(len(gene_order)):
        rv, fv = rho_matrix.values[i, j], fdr_matrix.values[i, j]
        if np.isnan(rv):
            continue
        star = "***" if fv < 0.001 else "**" if fv < 0.01 else "*" if fv < 0.05 else ""
        color = "white" if abs(rv) > vmax * 0.55 else "black"
        ax0.text(j, i, f"{rv:.2f}{star}", ha="center", va="center", fontsize=10, color=color)
cbar = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)
cbar.set_label("Spearman ρ (vs. MIF)", fontsize=10)
ax0.set_title("A. MIF co-expression across cancers", fontsize=12, loc="left")

# ---- Panel B: raw vs partial bars, with 95% CI error bars and bootstrap significance brackets ----
ax1 = axes[1]
x = np.arange(len(cancer_order)); width = 0.35

raw_vals = partial_ordered["raw_rho"].values
part_vals = partial_ordered["partial_rho"].values
ns = partial_ordered["n"].values

raw_ci = [fisher_ci(raw_vals[i], ns[i]) for i in range(len(cancer_order))]
part_ci = [fisher_ci(part_vals[i], ns[i]) for i in range(len(cancer_order))]
raw_err = [[raw_vals[i]-raw_ci[i][0] for i in range(len(cancer_order))],
           [raw_ci[i][1]-raw_vals[i] for i in range(len(cancer_order))]]
part_err = [[part_vals[i]-part_ci[i][0] for i in range(len(cancer_order))],
            [part_ci[i][1]-part_vals[i] for i in range(len(cancer_order))]]

ax1.bar(x - width/2, raw_vals, width, yerr=raw_err, capsize=4, label="Raw ρ",
        color="#C00000", error_kw={"elinewidth": 1.2, "ecolor": "black"})
ax1.bar(x + width/2, part_vals, width, yerr=part_err, capsize=4, label="Partial ρ\n(purity-adjusted)",
        color="#4472C4", error_kw={"elinewidth": 1.2, "ecolor": "black"})

# Significance brackets from the paired bootstrap test (bootstrap_diff_test.csv)
for i, c in enumerate(cancer_order):
    p_val = boot_ordered.loc[c, "bootstrap_p"]
    lower_of_pair = min(raw_ci[i][0], part_ci[i][0])
    bracket_y = lower_of_pair - 0.035
    x1, x2 = x[i] - width/2, x[i] + width/2
    ax1.plot([x1, x1, x2, x2], [bracket_y+0.012, bracket_y, bracket_y, bracket_y+0.012], color="black", linewidth=1)
    label = sig_label(p_val)
    ax1.text((x1+x2)/2, bracket_y-0.018, label, ha="center", va="top", fontsize=9,
              fontstyle="italic" if label == "n.s." else "normal")

ax1.set_xticks(x); ax1.set_xticklabels(cancer_order, fontsize=11)
ax1.axhline(0, color="black", linewidth=0.8)
ax1.set_ylabel("Spearman ρ (SP1 vs. MIF), 95% CI", fontsize=10)
ax1.set_title("B. SP1-MIF correlation: raw vs. purity-adjusted", fontsize=12, loc="left")
ax1.set_ylim(bottom=min(raw_ci[i][0] for i in range(len(cancer_order))) - 0.09, top=0.03)
ax1.legend(fontsize=9, loc="upper right")

plt.tight_layout()
plt.savefig("Phase3_MIF_correlation_figure_final.png", dpi=300, bbox_inches="tight")
plt.savefig("Phase3_MIF_correlation_figure_final.eps", format="eps", bbox_inches="tight")
print("Saved Phase3_MIF_correlation_figure_final.png and .eps")

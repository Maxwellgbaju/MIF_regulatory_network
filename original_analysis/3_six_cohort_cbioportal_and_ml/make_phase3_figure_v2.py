import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

corr = pd.read_csv("MIF_axis_correlations_6cancers.csv")
partial = pd.read_csv("SP1_MIF_partial_correlation.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
corr["cancer"] = corr["study"].map(study_to_code)
partial["cancer"] = partial["study"].map(study_to_code)
cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["SP1", "HIF1A", "RELA"]

def fisher_ci(rho, n):
    z = np.arctanh(rho)
    se = 1 / np.sqrt(n - 3)
    zcrit = 1.959963984540054
    return np.tanh(z - zcrit*se), np.tanh(z + zcrit*se)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.1, 1]})

# --- Panel A: heatmap (unchanged) ---
rho_matrix = corr.pivot(index="cancer", columns="gene", values="rho").reindex(index=cancer_order, columns=gene_order)
fdr_matrix = corr.pivot(index="cancer", columns="gene", values="fdr").reindex(index=cancer_order, columns=gene_order)
ax0 = axes[0]
vmax = np.nanmax(np.abs(rho_matrix.values))
im = ax0.imshow(rho_matrix.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax0.set_xticks(range(len(gene_order))); ax0.set_xticklabels(gene_order, fontsize=11)
ax0.set_yticks(range(len(cancer_order))); ax0.set_yticklabels(cancer_order, fontsize=11)
for i in range(len(cancer_order)):
    for j in range(len(gene_order)):
        rho_val = rho_matrix.values[i, j]; fdr_val = fdr_matrix.values[i, j]
        if np.isnan(rho_val): continue
        star = "***" if fdr_val < 0.001 else "**" if fdr_val < 0.01 else "*" if fdr_val < 0.05 else ""
        color = "white" if abs(rho_val) > vmax*0.55 else "black"
        ax0.text(j, i, f"{rho_val:.2f}{star}", ha="center", va="center", fontsize=10, color=color)
cbar = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)
cbar.set_label("Spearman ρ (vs. MIF)", fontsize=10)
ax0.set_title("A. MIF co-expression across cancers", fontsize=12, loc="left")

# --- Panel B: raw vs partial rho WITH 95% CI error bars ---
partial_ordered = partial.set_index("cancer").reindex(cancer_order).reset_index()
raw_lo, raw_hi = zip(*partial_ordered.apply(lambda r: fisher_ci(r["raw_rho"], r["n"]), axis=1))
part_lo, part_hi = zip(*partial_ordered.apply(lambda r: fisher_ci(r["partial_rho"], r["n"]), axis=1))
raw_err = [partial_ordered["raw_rho"] - np.array(raw_lo), np.array(raw_hi) - partial_ordered["raw_rho"]]
part_err = [partial_ordered["partial_rho"] - np.array(part_lo), np.array(part_hi) - partial_ordered["partial_rho"]]

ax1 = axes[1]
x = np.arange(len(cancer_order)); width = 0.35
ax1.bar(x - width/2, partial_ordered["raw_rho"], width, yerr=raw_err, capsize=4,
        label="Raw ρ", color="#C00000", error_kw={"elinewidth":1.2, "ecolor":"black"})
ax1.bar(x + width/2, partial_ordered["partial_rho"], width, yerr=part_err, capsize=4,
        label="Partial ρ\n(purity-adjusted)", color="#4472C4", error_kw={"elinewidth":1.2, "ecolor":"black"})
ax1.set_xticks(x); ax1.set_xticklabels(cancer_order, fontsize=11)
ax1.axhline(0, color="black", linewidth=0.8)
ax1.set_ylabel("Spearman ρ (SP1 vs. MIF), 95% CI", fontsize=10)
ax1.set_title("B. SP1-MIF correlation: raw vs. purity-adjusted", fontsize=12, loc="left")
ax1.legend(fontsize=9, loc="lower right")

plt.tight_layout()
plt.savefig("Phase3_MIF_correlation_figure_v2.png", dpi=300, bbox_inches="tight")
plt.savefig("Phase3_MIF_correlation_figure_v2.eps", format="eps", bbox_inches="tight")
print("Saved Phase3_MIF_correlation_figure_v2.png and .eps")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ---- Load data ----
corr = pd.read_csv("MIF_axis_correlations_6cancers.csv")
partial = pd.read_csv("SP1_MIF_partial_correlation.csv")

# Map long study IDs to short cancer codes for readable axis labels
study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC",
    "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD",
    "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC",
    "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
corr["cancer"] = corr["study"].map(study_to_code)
partial["cancer"] = partial["study"].map(study_to_code)

cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["SP1", "HIF1A", "RELA"]

# ---- Panel A: heatmap of rho, genes x cancers, FDR-significant cells starred ----
rho_matrix = corr.pivot(index="cancer", columns="gene", values="rho").reindex(index=cancer_order, columns=gene_order)
fdr_matrix = corr.pivot(index="cancer", columns="gene", values="fdr").reindex(index=cancer_order, columns=gene_order)

fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), gridspec_kw={"width_ratios": [1.1, 1]})

ax0 = axes[0]
vmax = np.nanmax(np.abs(rho_matrix.values))
im = ax0.imshow(rho_matrix.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

ax0.set_xticks(range(len(gene_order)))
ax0.set_xticklabels(gene_order, fontsize=11)
ax0.set_yticks(range(len(cancer_order)))
ax0.set_yticklabels(cancer_order, fontsize=11)

for i in range(len(cancer_order)):
    for j in range(len(gene_order)):
        rho_val = rho_matrix.values[i, j]
        fdr_val = fdr_matrix.values[i, j]
        if np.isnan(rho_val):
            continue
        star = ""
        if fdr_val < 0.001:
            star = "***"
        elif fdr_val < 0.01:
            star = "**"
        elif fdr_val < 0.05:
            star = "*"
        text_color = "white" if abs(rho_val) > vmax * 0.55 else "black"
        ax0.text(j, i, f"{rho_val:.2f}{star}", ha="center", va="center",
                  fontsize=10, color=text_color)

cbar = fig.colorbar(im, ax=ax0, fraction=0.046, pad=0.04)
cbar.set_label("Spearman ρ (vs. MIF)", fontsize=10)
ax0.set_title("A. MIF co-expression across cancers", fontsize=12, loc="left")

# ---- Panel B: raw vs partial rho for SP1, per cancer ----
ax1 = axes[1]
partial_ordered = partial.set_index("cancer").reindex(cancer_order).reset_index()

x = np.arange(len(cancer_order))
width = 0.35
bars1 = ax1.bar(x - width/2, partial_ordered["raw_rho"], width, label="Raw ρ", color="#C00000")
bars2 = ax1.bar(x + width/2, partial_ordered["partial_rho"], width, label="Partial ρ\n(purity-adjusted)", color="#4472C4")

ax1.set_xticks(x)
ax1.set_xticklabels(cancer_order, fontsize=11)
ax1.axhline(0, color="black", linewidth=0.8)
ax1.set_ylabel("Spearman ρ (SP1 vs. MIF)", fontsize=10)
ax1.set_title("B. SP1-MIF correlation: raw vs. purity-adjusted", fontsize=12, loc="left")
ax1.legend(fontsize=9, loc="lower right")

plt.tight_layout()
plt.savefig("Phase3_MIF_correlation_figure.png", dpi=300, bbox_inches="tight")
plt.savefig("Phase3_MIF_correlation_figure.eps", format="eps", bbox_inches="tight")
print("Saved Phase3_MIF_correlation_figure.png and .eps")

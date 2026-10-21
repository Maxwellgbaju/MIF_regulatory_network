import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = pd.read_csv("MIF_axis_immune_correlations_partial.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
res["cancer"] = res["study"].map(study_to_code)
cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["MIF", "SP1", "HIF1A", "RELA"]
markers = ["PDCD1", "FOXP3", "CD163"]

# ---- Panel A: three heatmaps (one per immune marker) of PARTIAL correlations ----
fig = plt.figure(figsize=(18, 6))
gs = fig.add_gridspec(1, 4, width_ratios=[1, 1, 1, 0.08], wspace=0.35)

vmax = np.nanmax(np.abs(res["partial_rho_ctrl_PTPRC"]))

for i, marker in enumerate(markers):
    ax = fig.add_subplot(gs[0, i])
    sub = res[res["immune_marker"] == marker]
    mat = sub.pivot(index="cancer", columns="axis_gene", values="partial_rho_ctrl_PTPRC").reindex(index=cancer_order, columns=gene_order)
    fdr_mat = sub.pivot(index="cancer", columns="axis_gene", values="partial_fdr").reindex(index=cancer_order, columns=gene_order)
    im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
    ax.set_xticks(range(len(gene_order))); ax.set_xticklabels(gene_order, fontsize=10)
    ax.set_yticks(range(len(cancer_order)))
    ax.set_yticklabels(cancer_order if i == 0 else [], fontsize=10)
    for r in range(len(cancer_order)):
        for c in range(len(gene_order)):
            v, f = mat.values[r, c], fdr_mat.values[r, c]
            if np.isnan(v): continue
            star = "***" if f < 0.001 else "**" if f < 0.01 else "*" if f < 0.05 else ""
            color = "white" if abs(v) > vmax * 0.55 else "black"
            ax.text(c, r, f"{v:.2f}{star}", ha="center", va="center", fontsize=8, color=color)
    ax.set_title(f"{marker}", fontsize=12)

cax = fig.add_subplot(gs[0, 3])
cbar = fig.colorbar(im, cax=cax)
cbar.set_label("Partial Spearman ρ\n(controlling for PTPRC)", fontsize=9)

plt.suptitle("A. MIF-axis gene correlation with immune markers, adjusted for overall immune infiltration", fontsize=13, y=1.03)
plt.savefig("Immune_partial_heatmaps.png", dpi=300, bbox_inches="tight")
plt.savefig("Immune_partial_heatmaps.eps", format="eps", bbox_inches="tight")
print("Saved Immune_partial_heatmaps.png/.eps")
plt.close()

# ---- Panel B: slope plot for PDCD1 specifically, showing the MIF/RELA vs SP1/HIF1A split ----
pdcd1 = res[res["immune_marker"] == "PDCD1"].copy()
pdcd1["cancer"] = pd.Categorical(pdcd1["cancer"], categories=cancer_order, ordered=True)
pdcd1 = pdcd1.sort_values(["axis_gene", "cancer"])

gene_colors = {"MIF": "#C00000", "RELA": "#ED7D31", "SP1": "#4472C4", "HIF1A": "#7030A0"}

fig, axes = plt.subplots(1, 4, figsize=(16, 5), sharey=True)
for ax, gene in zip(axes, gene_order):
    sub = pdcd1[pdcd1["axis_gene"] == gene]
    for _, row in sub.iterrows():
        ax.plot([0, 1], [row["raw_rho"], row["partial_rho_ctrl_PTPRC"]], marker="o",
                color=gene_colors[gene], alpha=0.8, linewidth=1.5)
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Raw", "Partial\n(ctrl. PTPRC)"], fontsize=10)
    ax.set_title(f"{gene} vs PDCD1", fontsize=12, color=gene_colors[gene])
    ax.set_xlim(-0.2, 1.2)

axes[0].set_ylabel("Spearman ρ", fontsize=11)
plt.suptitle("B. Raw-to-partial shift in correlation with PDCD1 (exhausted CD8+ T-cell marker), by cancer", fontsize=13)
plt.tight_layout()
plt.savefig("PDCD1_slopeplot.png", dpi=300, bbox_inches="tight")
plt.savefig("PDCD1_slopeplot.eps", format="eps", bbox_inches="tight")
print("Saved PDCD1_slopeplot.png/.eps")

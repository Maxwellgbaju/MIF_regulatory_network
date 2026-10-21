import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

res = pd.read_csv("MIF_axis_immune_correlations_partial.csv")
boot = pd.read_csv("immune_bootstrap_diff_all.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
res["cancer"] = res["study"].map(study_to_code)
boot["cancer"] = boot["study"].map(study_to_code)

cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["MIF", "SP1", "HIF1A", "RELA"]
cancer_colors = dict(zip(cancer_order, ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]))

pdcd1 = res[res["immune_marker"] == "PDCD1"].merge(
    boot[boot["immune_marker"] == "PDCD1"][["study", "axis_gene", "bootstrap_p"]],
    on=["study", "axis_gene"], how="left"
)

fig, axes = plt.subplots(1, 4, figsize=(17, 5.5), sharey=True)
for ax, gene in zip(axes, gene_order):
    sub = pdcd1[pdcd1["axis_gene"] == gene]
    for cancer in cancer_order:
        row = sub[sub["cancer"] == cancer]
        if row.empty:
            continue
        row = row.iloc[0]
        sig = row["bootstrap_p"] < 0.05
        ax.plot([0, 1], [row["raw_rho"], row["partial_rho_ctrl_PTPRC"]],
                marker="o", color=cancer_colors[cancer],
                linewidth=2.2 if sig else 1.0,
                linestyle="-" if sig else "--",
                alpha=0.95 if sig else 0.5)
    ax.axhline(0, color="black", linewidth=0.6, linestyle=":")
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Raw", "Partial\n(ctrl. PTPRC)"], fontsize=10)
    ax.set_title(f"{gene} vs PDCD1", fontsize=12)
    ax.set_xlim(-0.2, 1.2)

axes[0].set_ylabel("Spearman ρ", fontsize=11)

legend_elements = [Line2D([0], [0], color=cancer_colors[c], lw=2, label=c) for c in cancer_order]
legend_elements += [
    Line2D([0], [0], color="gray", lw=2.2, linestyle="-", label="Bootstrap p < 0.05 (shift significant)"),
    Line2D([0], [0], color="gray", lw=1.0, linestyle="--", alpha=0.6, label="Not significant"),
]
fig.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, -0.13), ncol=8, fontsize=9, frameon=False)

plt.suptitle("Raw-to-partial shift in correlation with PDCD1 (exhausted CD8+ T-cell marker), by cancer", fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig("PDCD1_slopeplot_v2.png", dpi=300, bbox_inches="tight")
plt.savefig("PDCD1_slopeplot_v2.eps", format="eps", bbox_inches="tight")
print("Saved PDCD1_slopeplot_v2.png/.eps")

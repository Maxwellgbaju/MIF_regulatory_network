import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = pd.read_csv("MIF_axis_survival_results_FDR.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
res["cancer"] = res["study"].map(study_to_code)
res["label"] = res["cancer"] + " - " + res["gene"]

# Order: group by cancer, then gene, matching your table order
cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["MIF", "SP1", "HIF1A", "RELA"]
res["cancer"] = pd.Categorical(res["cancer"], categories=cancer_order, ordered=True)
res["gene"] = pd.Categorical(res["gene"], categories=gene_order, ordered=True)
res = res.sort_values(["cancer", "gene"]).reset_index(drop=True)
res = res.iloc[::-1].reset_index(drop=True)  # reverse so KIRC appears at top of plot

y = np.arange(len(res))
significant = res["cox_fdr"] < 0.05

fig, ax = plt.subplots(figsize=(8, 9))

for i in range(len(res)):
    color = "#C00000" if significant.iloc[i] else "#888888"
    ax.plot([res["cox_HR_ci_lo"].iloc[i], res["cox_HR_ci_hi"].iloc[i]], [y[i], y[i]], color=color, linewidth=1.5, zorder=1)
    ax.scatter(res["cox_HR"].iloc[i], y[i], color=color, s=50 if significant.iloc[i] else 30, zorder=2,
               marker="D" if significant.iloc[i] else "o")

ax.axvline(1, color="black", linestyle="--", linewidth=1, zorder=0)
ax.set_yticks(y)
ax.set_yticklabels(res["label"], fontsize=9)
ax.set_xlabel("Hazard ratio (95% CI), continuous expression Cox model", fontsize=10)
ax.set_xscale("log")
ax.set_title("A. Overall survival association across 24 gene-cancer tests", fontsize=12, loc="left")

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker="D", color="w", markerfacecolor="#C00000", markersize=8, label="FDR < 0.05"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="#888888", markersize=7, label="Not significant"),
]
ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

plt.tight_layout()
plt.savefig("Forest_plot_survival_24tests.png", dpi=300, bbox_inches="tight")
plt.savefig("Forest_plot_survival_24tests.eps", format="eps", bbox_inches="tight")
print("Saved Forest_plot_survival_24tests.png and .eps")

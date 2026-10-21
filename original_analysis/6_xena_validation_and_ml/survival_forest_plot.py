import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("full_survival_screen_results.csv")
sig = df[df["FDR"] < 0.05].copy()
sig["label"] = sig["gene"] + " — " + sig["cancer_type"]
sig = sig.sort_values("cox_HR")

# Flag the Brain LGG cluster explicitly (5 genes, same cancer type, same direction --
# a potential shared-confound signature, not necessarily 5 independent findings)
sig["is_brain_lgg"] = sig["cancer_type"] == "Brain Lower Grade Glioma"

fig, ax = plt.subplots(figsize=(10, 9))
y = np.arange(len(sig))

colors = []
for _, row in sig.iterrows():
    if row["is_brain_lgg"]:
        colors.append("#7030A0")  # distinct color to flag the Brain LGG cluster
    elif row["cox_HR"] > 1:
        colors.append("#C00000")  # harmful
    else:
        colors.append("#2E7D32")  # protective

for i, (_, row) in enumerate(sig.iterrows()):
    ax.plot([row["cox_HR_lo"], row["cox_HR_hi"]], [i, i], color=colors[i], linewidth=1.8)
    marker = "D" if row["is_brain_lgg"] else "o"
    ax.scatter(row["cox_HR"], i, color=colors[i], s=70, zorder=3, marker=marker)

ax.axvline(1, color="black", linestyle="--", linewidth=1)
ax.set_yticks(y)
ax.set_yticklabels(sig["label"], fontsize=9)
ax.set_xscale("log")
ax.set_xlabel("Hazard Ratio (95% CI), log scale", fontsize=10)
ax.set_title("Significant survival associations (FDR < 0.05)\n" +
             "UCSC Xena TCGA TARGET GTEx cohort, 403 gene-by-cancer-type tests\n" +
             "Diamonds = Brain Lower Grade Glioma cluster (possible shared confound, see caption)",
             fontsize=11)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Patch(facecolor="#C00000", label="Harmful (HR > 1)"),
    Patch(facecolor="#2E7D32", label="Protective (HR < 1)"),
    Line2D([0], [0], marker="D", color="w", markerfacecolor="#7030A0", markersize=9, label="Brain LGG cluster"),
]
ax.legend(handles=legend_elements, loc="lower right", fontsize=9)

plt.tight_layout()
plt.savefig("survival_forest_plot.png", dpi=300, bbox_inches="tight")
plt.savefig("survival_forest_plot.eps", format="eps", bbox_inches="tight")
print("Saved survival_forest_plot.png/.eps")

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

cbio_stab = pd.read_csv("cbioportal_stability_selection.csv")
xena_stab = pd.read_csv("/shared/rc/malariarna/MIF_student_network_v2/gepia_test/xena_stability_selection.csv")
rsf = pd.read_csv("capstone_rsf_importance.csv")

cbio_stab = cbio_stab.set_index("gene")["selection_pct"]
xena_stab = xena_stab.set_index("gene")["selection_pct"]

all_genes = ["MIF", "RELA", "EGR1", "SP1", "NR1H2", "RXRA", "IRF1", "KLF4", "IRF2",
             "NFIL3", "ESR1", "TP53", "TFAP2A", "NORAD", "AL356488.2"]

convergent = {"TFAP2A", "KLF4", "NORAD", "IRF1", "TP53"}

fig, axes = plt.subplots(1, 3, figsize=(19, 7.5))

def get_color(g):
    if g in convergent:
        return "#C00000"
    if g == "RELA":
        return "#4472C4"
    if g == "RXRA":
        return "#E8A33D"
    return "#BBBBBB"

genes_a = [g for g in all_genes if g in cbio_stab.index]
vals_a = [cbio_stab[g] for g in genes_a]
order = np.argsort(vals_a)
genes_a = [genes_a[i] for i in order]
vals_a = [vals_a[i] for i in order]
axes[0].barh(genes_a, vals_a, color=[get_color(g) for g in genes_a], edgecolor="black", linewidth=0.6)
axes[0].set_title("A. cBioPortal LASSO\nStability selection (n=510, 100 resampled fits)", fontsize=11)
axes[0].set_xlabel("% of fits selecting this gene", fontsize=10)
axes[0].axvline(80, color="black", linestyle="--", linewidth=0.8, alpha=0.5)

genes_b = [g for g in all_genes if g in xena_stab.index]
vals_b = [xena_stab[g] for g in genes_b]
order = np.argsort(vals_b)
genes_b = [genes_b[i] for i in order]
vals_b = [vals_b[i] for i in order]
axes[1].barh(genes_b, vals_b, color=[get_color(g) for g in genes_b], edgecolor="black", linewidth=0.6)
axes[1].set_title("B. Xena LASSO\nStability selection (n=530, 100 resampled fits)", fontsize=11)
axes[1].set_xlabel("% of fits selecting this gene", fontsize=10)
axes[1].axvline(80, color="black", linestyle="--", linewidth=0.8, alpha=0.5)

genes_c = [g for g in all_genes if g in rsf["gene"].values]
rsf_idx = rsf.set_index("gene")
vals_c = [rsf_idx.loc[g, "importance_mean"] for g in genes_c]
errs_c = [rsf_idx.loc[g, "importance_std"] for g in genes_c]
order = np.argsort(vals_c)
genes_c = [genes_c[i] for i in order]
vals_c = [vals_c[i] for i in order]
errs_c = [errs_c[i] for i in order]
axes[2].barh(genes_c, vals_c, xerr=errs_c, color=[get_color(g) for g in genes_c], edgecolor="black", linewidth=0.6, capsize=3)
axes[2].set_title("C. Random Survival Forest\nPermutation importance +/- SD (15 repeats)", fontsize=11)
axes[2].set_xlabel("Importance (mean C-index drop)", fontsize=10)

legend_elements = [
    Patch(facecolor="#C00000", label="Cross-database convergent"),
    Patch(facecolor="#4472C4", label="RELA (cBioPortal-specific)"),
    Patch(facecolor="#E8A33D", label="RXRA (Xena-specific, unexplained)"),
    Patch(facecolor="#BBBBBB", label="Other"),
]
fig.legend(handles=legend_elements, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.05), fontsize=10)
plt.suptitle("Capstone: Cross-database, cross-method machine learning risk modeling (KIRC)", fontsize=13, y=1.02)
plt.tight_layout()
plt.savefig("capstone_final_figure.png", dpi=300, bbox_inches="tight")
plt.savefig("capstone_final_figure.eps", format="eps", bbox_inches="tight")
print("Saved capstone_final_figure.png/.eps")

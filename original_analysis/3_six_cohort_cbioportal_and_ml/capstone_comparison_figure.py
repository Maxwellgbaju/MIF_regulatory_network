import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

cbio = pd.read_csv("capstone_cbioportal_lasso_coefficients.csv").set_index("gene")
xena = pd.read_csv("/shared/rc/malariarna/MIF_student_network_v2/gepia_test/capstone_lasso_coefficients.csv").set_index("gene")
rsf = pd.read_csv("capstone_rsf_importance.csv").set_index("gene")

all_genes = ["MIF", "RELA", "EGR1", "SP1", "NR1H2", "RXRA", "IRF1", "KLF4", "IRF2",
             "NFIL3", "ESR1", "TP53", "TFAP2A", "NORAD", "AL356488.2"]

cross_db_convergent = {"IRF1", "KLF4", "TP53", "TFAP2A", "NORAD"}

fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharey=False)

def plot_panel(ax, values_dict, title, ylabel):
    genes_here = [g for g in all_genes if g in values_dict.index]
    vals = [values_dict.loc[g].iloc[0] if hasattr(values_dict.loc[g], 'iloc') else values_dict.loc[g] for g in genes_here]
    colors = ["#C00000" if g in cross_db_convergent else ("#4472C4" if g == "RELA" else "#BBBBBB") for g in genes_here]
    order = np.argsort(vals)
    genes_sorted = [genes_here[i] for i in order]
    vals_sorted = [vals[i] for i in order]
    colors_sorted = [colors[i] for i in order]
    ax.barh(genes_sorted, vals_sorted, color=colors_sorted, edgecolor="black", linewidth=0.6)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_title(title, fontsize=12)
    ax.set_xlabel(ylabel, fontsize=10)

cbio_col = cbio["standardized_coefficient"]
xena_col = xena["standardized_coefficient"]
rsf_col = rsf["importance_mean"] if "importance_mean" in rsf.columns else rsf["importance"]

plot_panel(axes[0], cbio_col.to_frame(), "A. cBioPortal LASSO\n(n=510)", "Standardized coefficient")
plot_panel(axes[1], xena_col.to_frame(), "B. Xena LASSO\n(n=530)", "Standardized coefficient")
plot_panel(axes[2], rsf_col.to_frame(), "C. Random Survival Forest\n(cBioPortal, permutation importance)", "Importance")

from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#C00000", label="Cross-database convergent (5 genes)"),
    Patch(facecolor="#4472C4", label="RELA"),
    Patch(facecolor="#BBBBBB", label="Other"),
]
fig.legend(handles=legend_elements, loc="lower center", ncol=3, bbox_to_anchor=(0.5, -0.05), fontsize=10)

plt.suptitle("Capstone: Machine learning risk modeling across two independently processed cohorts\n(KIRC, MIF-axis regulatory candidates)", fontsize=13, y=1.03)
plt.tight_layout()
plt.savefig("capstone_ml_comparison.png", dpi=300, bbox_inches="tight")
plt.savefig("capstone_ml_comparison.eps", format="eps", bbox_inches="tight")
print("Saved capstone_ml_comparison.png/.eps")

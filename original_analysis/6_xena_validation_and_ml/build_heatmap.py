import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("tumor_normal_full_results.csv")

# Rank cancer types by MIF's own FDR (the gene this whole study is about)
mif_ranked = df[df["gene"] == "MIF"].sort_values("FDR")
top10_cancers = mif_ranked.head(10)["cancer_type"].tolist()

# Explicitly add Pancreatic Adenocarcinoma for clinical severity (deadliest common solid tumor),
# even though it ranks #11 by MIF significance -- not silently, but as a stated addition
clinical_addition = "Pancreatic Adenocarcinoma"
selected_cancers = top10_cancers + [clinical_addition]
print(f"Top 10 by MIF significance: {top10_cancers}")
print(f"Clinical addition: {clinical_addition}")

genes = ["MIF", "SP1", "TP53", "EGR1", "KLF4", "IRF1", "IRF2", "NFIL3", "ESR1", "TFAP2A", "NR1H2", "RXRA"]

sub = df[df["cancer_type"].isin(selected_cancers)]
fc_matrix = sub.pivot(index="gene", columns="cancer_type", values="log2FC").reindex(index=genes, columns=selected_cancers)
fdr_matrix = sub.pivot(index="gene", columns="cancer_type", values="FDR").reindex(index=genes, columns=selected_cancers)

def sig_label(fdr):
    if pd.isna(fdr): return ""
    if fdr < 0.001: return "***"
    if fdr < 0.01: return "**"
    if fdr < 0.05: return "*"
    return "ns"

fig, ax = plt.subplots(figsize=(14, 8))
vmax = np.nanmax(np.abs(fc_matrix.values))
im = ax.imshow(fc_matrix.values, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

ax.set_xticks(range(len(selected_cancers)))
short_labels = [c if c != clinical_addition else f"{c}\n(clinical addition)" for c in selected_cancers]
ax.set_xticklabels(short_labels, rotation=40, ha="right", fontsize=9)
ax.set_yticks(range(len(genes)))
ax.set_yticklabels(genes, fontsize=10, fontweight="bold")

for i in range(len(genes)):
    for j in range(len(selected_cancers)):
        fc = fc_matrix.values[i, j]
        fdr = fdr_matrix.values[i, j]
        if pd.isna(fc):
            continue
        label = sig_label(fdr)
        color = "white" if abs(fc) > vmax * 0.55 else "black"
        ax.text(j, i, f"{fc:.2f}\n{label}", ha="center", va="center", fontsize=7.5, color=color)

# Mark the clinical-addition column with a box
clinical_idx = selected_cancers.index(clinical_addition)
ax.axvline(clinical_idx - 0.5, color="black", linewidth=2)
ax.axvline(clinical_idx + 0.5, color="black", linewidth=2)

cbar = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
cbar.set_label("log2 fold-change (Tumor vs. Normal)", fontsize=10)

ax.set_title("Pan-cancer tumor-vs-normal expression: MIF and 11 candidate transcription factors\n" +
             "Top 10 cancer types by MIF significance, plus Pancreatic Adenocarcinoma (clinical addition, outlined)\n" +
             "*** FDR<0.001, ** FDR<0.01, * FDR<0.05, ns = not significant", fontsize=11)

plt.tight_layout()
plt.savefig("MIF_pancancer_heatmap.png", dpi=300, bbox_inches="tight")
plt.savefig("MIF_pancancer_heatmap.eps", format="eps", bbox_inches="tight")
print("Saved MIF_pancancer_heatmap.png/.eps")

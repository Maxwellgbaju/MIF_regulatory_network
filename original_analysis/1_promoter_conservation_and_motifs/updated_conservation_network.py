import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(13, 13))
ax.set_xlim(-7, 7)
ax.set_ylim(-7, 7)
ax.axis("off")

# MIF center
ax.scatter([0], [0], s=3200, color="#D32F2F", edgecolor="black", linewidth=2, zorder=5)
ax.text(0, 0, "MIF", ha="center", va="center", fontsize=15, fontweight="bold", color="white", zorder=6)

# Conservation tier colors
tier_colors = {
    "full3": "#2E7D32",       # dark green: full 3-species conservation
    "weak3": "#8BC34A",       # light green: present all 3, sub-threshold
    "heterodimer": "#FFA726", # orange: heterodimer-only representation
    "single": "#FFEB3B",      # yellow: single-species only
    "none": "#BDBDBD",        # grey: below threshold everywhere
}

genes = {
    "SP1":    ("full3", True),   "HIF1A": ("full3", True),  "RELA": ("full3", True),
    "EGR1":   ("full3", False),  "KLF4": ("full3", False),  "NFIL3": ("full3", False), "TFAP2A": ("full3", False),
    "NR1H2":  ("weak3", False),
    "RXRA":   ("heterodimer", False),
    "IRF1":   ("single", False),
    "TP53":   ("none", False), "ESR1": ("none", False), "IRF2": ("none", False),
}

n = len(genes)
for i, (gene, (tier, curated)) in enumerate(genes.items()):
    angle = 2 * np.pi * i / n
    x, y = 4.8 * np.cos(angle), 4.8 * np.sin(angle)
    color = tier_colors[tier]
    edge_width = 3.5 if curated else 1.2
    edge_style = "-"
    ax.scatter([x], [y], s=1700, color=color, edgecolor="black", linewidth=edge_width, zorder=4)
    ax.text(x, y, gene, ha="center", va="center", fontsize=9.5, fontweight="bold" if curated else "normal", zorder=5)
    ax.plot([0, x], [0, y], color="#888888", linewidth=1, alpha=0.5, zorder=1)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D
legend_elements = [
    Patch(facecolor=tier_colors["full3"], edgecolor="black", label="Full 3-species conservation"),
    Patch(facecolor=tier_colors["weak3"], edgecolor="black", label="Present all 3 species, sub-threshold"),
    Patch(facecolor=tier_colors["heterodimer"], edgecolor="black", label="Heterodimer motif only (RXRA)"),
    Patch(facecolor=tier_colors["single"], edgecolor="black", label="Single-species match only"),
    Patch(facecolor=tier_colors["none"], edgecolor="black", label="Below threshold, all 3 species"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="white", markeredgecolor="black",
           markeredgewidth=3.5, markersize=14, label="Thick border = original curated panel (SP1/HIF1A/RELA)"),
]
ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(-0.15, 1.12), fontsize=9, frameon=True)
ax.set_title("MIF candidate transcription factors: conservation tier\n(14 genes -- original curated panel + broader JASPAR-score screen)",
              fontsize=13, pad=20)

plt.tight_layout()
plt.savefig("MIF_TF_conservation_network.png", dpi=300, bbox_inches="tight")
plt.savefig("MIF_TF_conservation_network.eps", format="eps", bbox_inches="tight")
print("Saved MIF_TF_conservation_network.png/.eps")

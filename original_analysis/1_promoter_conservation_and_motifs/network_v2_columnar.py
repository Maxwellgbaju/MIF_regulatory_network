import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import numpy as np

fig, ax = plt.subplots(figsize=(18, 13))
ax.set_xlim(-2, 20)
ax.set_ylim(-1, 14)
ax.axis("off")

def curve(ax, p0, p1, color, lw, style="-", curvature=0.3):
    x0, y0 = p0; x1, y1 = p1
    mx, my = (x0+x1)/2, (y0+y1)/2
    dx, dy = x1-x0, y1-y0
    nx, ny = -dy, dx
    norm = np.hypot(nx, ny) + 1e-9
    ctrl = (mx + nx/norm*curvature*np.hypot(dx,dy), my + ny/norm*curvature*np.hypot(dx,dy))
    path = Path([p0, ctrl, p1], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
    patch = mpatches.PathPatch(path, facecolor="none", edgecolor=color, lw=lw, linestyle=style, zorder=1)
    ax.add_patch(patch)

# ---- Column x-positions ----
X_TF, X_MIF, X_MI, X_LNC = 0, 6, 12, 18

# ---- TF column (14 genes, tier as text) ----
tfs = [
    ("SP1", "4-tier"), ("RELA", "3-tier"), ("HIF1A", "3-tier"),
    ("TFAP2A", "3-species"), ("EGR1", "3-species"), ("KLF4", "3-species"), ("NFIL3", "3-species"),
    ("NR1H2", "sub-threshold"), ("RXRA", "heterodimer only"), ("IRF1", "mouse only"),
    ("TP53", "below threshold"), ("ESR1", "below threshold"), ("IRF2", "below threshold"),
]
n_tf = len(tfs)
tf_y = np.linspace(13, 1, n_tf)
tf_pos = {}
for (gene, tier), y in zip(tfs, tf_y):
    tf_pos[gene] = (X_TF, y)
    ax.scatter([X_TF], [y], s=1600, color="#4472C4", edgecolor="black", linewidth=1.3, zorder=3)
    ax.text(X_TF, y, gene, ha="center", va="center", fontsize=8, fontweight="bold", color="white", zorder=4)
    ax.text(X_TF, y-0.55, f"({tier})", ha="center", va="center", fontsize=6.8, color="#333333", zorder=4)

# ---- MIF ----
mif_y = 7
ax.scatter([X_MIF], [mif_y], s=3200, color="#C00000", edgecolor="black", linewidth=1.5, zorder=3)
ax.text(X_MIF, mif_y, "MIF", ha="center", va="center", fontsize=13, fontweight="bold", color="white", zorder=4)

for gene, (x, y) in tf_pos.items():
    curve(ax, (x+0.35, y), (X_MIF-0.55, mif_y), "#4472C4", 1.4, "-", curvature=0.15)

# ---- miRNA column ----
mirnas_original = ["hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p", "hsa-miR-23a-5p",
                    "hsa-miR-363-5p", "hsa-miR-23b-5p", "hsa-miR-541-3p", "hsa-miR-92a-1-5p",
                    "hsa-miR-608", "hsa-miR-1537-3p"]
rediscovered = {"hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p"}
n_mi = len(mirnas_original)
mi_y = np.linspace(13, 1, n_mi)
mi_pos = {}
for mirna, y in zip(mirnas_original, mi_y):
    mi_pos[mirna] = (X_MI, y)
    color = "#70AD47" if mirna not in rediscovered else "#548235"
    ax.scatter([X_MI], [y], s=1300, color=color, edgecolor="black", linewidth=1.3, zorder=3)
    ax.text(X_MI, y, mirna, ha="center", va="center", fontsize=7.2, color="white", fontweight="bold", zorder=4)
    curve(ax, (X_MIF+0.55, mif_y), (X_MI-0.5, y), "#C00000", 1.3, "--", curvature=-0.15)

# ---- lncRNA column ----
lncrnas = {
    "AL356488.2": ["hsa-miR-654-5p"],
    "AC020922.3": ["hsa-miR-654-5p", "hsa-miR-541-3p"],
    "AL355987.4": ["hsa-miR-541-3p", "hsa-miR-23a-5p"],
    "LMNTD2-AS1": ["hsa-miR-541-3p", "hsa-miR-92a-1-5p"],
    "NORAD": ["hsa-miR-451a"],
}
n_lnc = len(lncrnas)
lnc_y = np.linspace(11, 3, n_lnc)
for (lnc, targets), y in zip(lncrnas.items(), lnc_y):
    is_flag = lnc == "AL356488.2"
    facecolor = "#FFC000" if not is_flag else "#FFD966"
    lw = 2.6 if is_flag else 1.3
    rect = mpatches.FancyBboxPatch((X_LNC-0.9, y-0.35), 1.8, 0.7, boxstyle="round,pad=0.03,rounding_size=0.06",
                                     facecolor=facecolor, edgecolor="black", linewidth=lw, zorder=3)
    ax.add_patch(rect)
    ax.text(X_LNC, y, lnc, ha="center", va="center", fontsize=7.5, fontweight="bold" if is_flag else "normal", zorder=4)
    for t in targets:
        if t in mi_pos:
            tx, ty = mi_pos[t]
            curve(ax, (X_LNC-0.9, y), (tx+0.45, ty), "#BF9000", 1.2, ":", curvature=0.1)

# ---- Legend ----
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], color="#4472C4", lw=2, label="TF -> MIF (activation)"),
    Line2D([0],[0], color="#C00000", lw=2, linestyle="--", label="miRNA -| MIF (repression)"),
    Line2D([0],[0], color="#BF9000", lw=2, linestyle=":", label="lncRNA -> miRNA (sponge)"),
    mpatches.Patch(facecolor="#548235", edgecolor="black", label="miRNA: independently rediscovered"),
    mpatches.Patch(facecolor="#FFD966", edgecolor="black", label="lncRNA: flagship (AL356488.2)"),
]
ax.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, -0.08), ncol=3, fontsize=9.5, frameon=True)
ax.set_title("MIF Multi-Layer Regulatory Network: Transcription Factors, MicroRNAs, and lncRNA Sponges", fontsize=14, pad=15)

plt.tight_layout()
plt.savefig("MIF_network_v2_columnar.png", dpi=300, bbox_inches="tight")
plt.savefig("MIF_network_v2_columnar.eps", format="eps", bbox_inches="tight")
print("Saved MIF_network_v2_columnar.png/.eps")

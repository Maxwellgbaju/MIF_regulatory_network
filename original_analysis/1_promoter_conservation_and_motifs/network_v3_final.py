import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import numpy as np

fig, ax = plt.subplots(figsize=(19, 13))
ax.set_xlim(-2, 20)
ax.set_ylim(-1, 14)
ax.axis("off")

def curve(ax, p0, p1, color, lw, style="-", curvature=0.3, alpha=1.0):
    x0, y0 = p0; x1, y1 = p1
    mx, my = (x0+x1)/2, (y0+y1)/2
    dx, dy = x1-x0, y1-y0
    nx, ny = -dy, dx
    norm = np.hypot(nx, ny) + 1e-9
    ctrl = (mx + nx/norm*curvature*np.hypot(dx,dy), my + ny/norm*curvature*np.hypot(dx,dy))
    path = Path([p0, ctrl, p1], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
    patch = mpatches.PathPatch(path, facecolor="none", edgecolor=color, lw=lw, linestyle=style, alpha=alpha, zorder=1)
    ax.add_patch(patch)

X_TF, X_MIF, X_MI, X_LNC = 0, 6, 12, 18

tier_colors = {"full3": "#2E7D32", "weak3": "#8BC34A", "heterodimer": "#FFA726", "single": "#FFEE58", "none": "#BDBDBD"}

tfs = [
    ("SP1", "full3", "4-tier", True), ("RELA", "full3", "3-tier", True), ("HIF1A", "full3", "3-tier", True),
    ("TFAP2A", "full3", "3-species", False), ("EGR1", "full3", "3-species", False),
    ("KLF4", "full3", "3-species", False), ("NFIL3", "full3", "3-species", False),
    ("NR1H2", "weak3", "sub-threshold", False), ("RXRA", "heterodimer", "heterodimer only", False),
    ("IRF1", "single", "mouse only", False),
    ("TP53", "none", "below threshold", False), ("ESR1", "none", "below threshold", False), ("IRF2", "none", "below threshold", False),
]
n_tf = len(tfs)
tf_y = np.linspace(13, 1, n_tf)
tf_pos = {}
for (gene, tier, label, curated), y in zip(tfs, tf_y):
    tf_pos[gene] = (X_TF, y)
    color = tier_colors[tier]
    lw = 3.0 if curated else 1.2
    ax.scatter([X_TF], [y], s=1600, color=color, edgecolor="black", linewidth=lw, zorder=3)
    ax.text(X_TF, y, gene, ha="center", va="center", fontsize=8, fontweight="bold", color="black", zorder=4)
    ax.text(X_TF, y-0.58, f"({label})", ha="center", va="center", fontsize=6.5, color="#333333", zorder=4)

mif_y = 7
ax.scatter([X_MIF], [mif_y], s=3200, color="#C00000", edgecolor="black", linewidth=1.5, zorder=3)
ax.text(X_MIF, mif_y, "MIF", ha="center", va="center", fontsize=13, fontweight="bold", color="white", zorder=4)

for gene, (x, y) in tf_pos.items():
    curve(ax, (x+0.35, y), (X_MIF-0.55, mif_y), "#4472C4", 1.4, "-", curvature=0.15)

mirnas_original = ["hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p", "hsa-miR-23a-5p",
                    "hsa-miR-363-5p", "hsa-miR-23b-5p", "hsa-miR-541-3p", "hsa-miR-92a-1-5p",
                    "hsa-miR-608", "hsa-miR-1537-3p"]
rediscovered = {"hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p"}
n_mi = len(mirnas_original)
mi_y = np.linspace(13, 1, n_mi)
mi_pos = {}
for mirna, y in zip(mirnas_original, mi_y):
    mi_pos[mirna] = (X_MI, y)
    color = "#A9D18E" if mirna not in rediscovered else "#548235"
    lw = 2.6 if mirna in rediscovered else 1.2
    ax.scatter([X_MI], [y], s=1300, color=color, edgecolor="black", linewidth=lw, zorder=3)
    ax.text(X_MI, y, mirna, ha="center", va="center", fontsize=7.2, color="black", fontweight="bold", zorder=4)
    curve(ax, (X_MIF+0.55, mif_y), (X_MI-0.5, y), "#C00000", 1.3, "--", curvature=-0.15)

# lncRNA -> miRNA edges; only hsa-miR-654-5p:AL356488.2 has a confirmed real CLIP-seq count (28 experiments)
lncrnas = {
    "AL356488.2": {"hsa-miR-654-5p": 28},
    "AC020922.3": {"hsa-miR-654-5p": None, "hsa-miR-541-3p": None},
    "AL355987.4": {"hsa-miR-541-3p": None, "hsa-miR-23a-5p": None},
    "LMNTD2-AS1": {"hsa-miR-541-3p": None, "hsa-miR-92a-1-5p": None},
    "NORAD": {"hsa-miR-451a": None},
}
n_lnc = len(lncrnas)
lnc_y = np.linspace(11, 3, n_lnc)
for (lnc, targets), y in zip(lncrnas.items(), lnc_y):
    is_flag = lnc == "AL356488.2"
    facecolor = "#FFD966" if not is_flag else "#FFC000"
    lw = 2.8 if is_flag else 1.3
    rect = mpatches.FancyBboxPatch((X_LNC-0.9, y-0.35), 1.8, 0.7, boxstyle="round,pad=0.03,rounding_size=0.06",
                                     facecolor=facecolor, edgecolor="black", linewidth=lw, zorder=3)
    ax.add_patch(rect)
    ax.text(X_LNC, y, lnc, ha="center", va="center", fontsize=7.5, fontweight="bold" if is_flag else "normal", zorder=4)
    for t, clip_n in targets.items():
        if t in mi_pos:
            tx, ty = mi_pos[t]
            edge_lw = 3.2 if clip_n else 1.4
            edge_color = "#BF6F00" if clip_n else "#D9A441"
            curve(ax, (X_LNC-0.9, y), (tx+0.45, ty), edge_color, edge_lw, "-", curvature=0.12, alpha=0.85)
            if clip_n:
                mx, my = (X_LNC-0.9+tx+0.45)/2, (y+ty)/2 + 0.15
                ax.text(mx, my, f"{clip_n} CLIP-seq exp.", fontsize=6.5, color="#8A5300", ha="center", zorder=5,
                        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.8))

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0],[0], color="#4472C4", lw=2, label="TF -> MIF (activation)"),
    Line2D([0],[0], color="#C00000", lw=2, linestyle="--", label="miRNA -| MIF (repression)"),
    Line2D([0],[0], color="#BF6F00", lw=3, label="lncRNA -> miRNA (CLIP-seq n given where confirmed)"),
    mpatches.Patch(facecolor=tier_colors["full3"], edgecolor="black", label="TF: full 3-species conserved"),
    mpatches.Patch(facecolor=tier_colors["weak3"], edgecolor="black", label="TF: sub-threshold, all species"),
    mpatches.Patch(facecolor=tier_colors["heterodimer"], edgecolor="black", label="TF: heterodimer motif only"),
    mpatches.Patch(facecolor=tier_colors["single"], edgecolor="black", label="TF: single-species only"),
    mpatches.Patch(facecolor=tier_colors["none"], edgecolor="black", label="TF: below threshold, all species"),
    mpatches.Patch(facecolor="#548235", edgecolor="black", label="miRNA: found by 2 independent tools"),
    mpatches.Patch(facecolor="#FFC000", edgecolor="black", label="lncRNA: flagship (AL356488.2)"),
]
ax.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, -0.11), ncol=3, fontsize=8.5, frameon=True)
ax.set_title("MIF Multi-Layer Regulatory Network: Transcription Factors, MicroRNAs, and lncRNA Sponges", fontsize=14, pad=15)
ax.text(X_TF, 13.7, "Thick border = original curated panel (SP1/RELA/HIF1A)", fontsize=8, ha="left", style="italic", color="#333333")

plt.tight_layout()
plt.savefig("MIF_network_v3_final.png", dpi=300, bbox_inches="tight")
plt.savefig("MIF_network_v3_final.eps", format="eps", bbox_inches="tight")
print("Saved MIF_network_v3_final.png/.eps")

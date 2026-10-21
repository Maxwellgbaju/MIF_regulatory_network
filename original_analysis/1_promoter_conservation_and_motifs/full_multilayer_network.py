import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig, ax = plt.subplots(figsize=(16, 14))
ax.set_xlim(-11, 11)
ax.set_ylim(-9, 9)
ax.axis("off")

# ---- MIF center ----
ax.scatter([0], [0], s=3600, color="#D32F2F", edgecolor="black", linewidth=2, zorder=6)
ax.text(0, 0, "MIF", ha="center", va="center", fontsize=16, fontweight="bold", color="white", zorder=7)

tier_colors = {
    "full3": "#2E7D32", "weak3": "#8BC34A", "heterodimer": "#FFA726",
    "single": "#FFEB3B", "none": "#BDBDBD",
}

# ---- TF layer: left semicircle, activation arrows pointing INTO MIF ----
tfs = {
    "SP1": ("full3", True), "HIF1A": ("full3", True), "RELA": ("full3", True),
    "EGR1": ("full3", False), "KLF4": ("full3", False), "NFIL3": ("full3", False), "TFAP2A": ("full3", False),
    "NR1H2": ("weak3", False), "RXRA": ("heterodimer", False), "IRF1": ("single", False),
    "TP53": ("none", False), "ESR1": ("none", False), "IRF2": ("none", False),
}
n_tf = len(tfs)
tf_pos = {}
for i, (gene, (tier, curated)) in enumerate(tfs.items()):
    angle = np.pi/2 + (np.pi * (i+0.5) / n_tf)  # spread across left half (90 to 270 deg)
    x, y = 7.2*np.cos(angle), 7.2*np.sin(angle)
    tf_pos[gene] = (x, y)
    color = tier_colors[tier]
    lw = 3.2 if curated else 1.1
    ax.scatter([x], [y], s=1500, color=color, edgecolor="black", linewidth=lw, zorder=4)
    ax.text(x, y, gene, ha="center", va="center", fontsize=8.5, fontweight="bold" if curated else "normal", zorder=5)
    # activation arrow: pointed arrowhead into MIF
    ax.annotate("", xy=(0,0), xytext=(x*0.82, y*0.82),
                arrowprops=dict(arrowstyle="-|>", color="#2E7D32", lw=1.4, mutation_scale=15), zorder=2)

# ---- miRNA layer: right semicircle, inhibition (T-bar) arrows into MIF ----
mirnas_original = ["hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p", "hsa-miR-23a-5p",
                    "hsa-miR-363-5p", "hsa-miR-23b-5p", "hsa-miR-541-3p", "hsa-miR-92a-1-5p",
                    "hsa-miR-608", "hsa-miR-1537-3p"]
rediscovered = {"hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p"}  # independently found again today

n_mi = len(mirnas_original)
mi_pos = {}
for i, mirna in enumerate(mirnas_original):
    angle = -np.pi/2 + (np.pi * (i+0.5) / n_mi)  # spread across right half (-90 to 90 deg)
    x, y = 7.2*np.cos(angle), 7.2*np.sin(angle)
    mi_pos[mirna] = (x, y)
    is_redis = mirna in rediscovered
    color = "#F5C242" if is_redis else "#FFE9A8"
    lw = 3.0 if is_redis else 1.0
    ax.scatter([x], [y], s=1150, color=color, edgecolor="black", linewidth=lw, zorder=4)
    ax.text(x, y, mirna, ha="center", va="center", fontsize=7, zorder=5)
    # inhibition arrow: T-bar into MIF
    dx, dy = x*0.85 - 0, y*0.85 - 0
    ax.plot([x*0.85, 0], [y*0.85, 0], color="#C00000", lw=1.4, zorder=2)
    # draw T-bar near MIF end
    ux, uy = -dx, -dy
    norm = np.hypot(ux, uy)
    ux, uy = ux/norm, uy/norm
    px, py = -uy, ux
    bx, by = ux*1.1, uy*1.1
    ax.plot([bx-px*0.35, bx+px*0.35], [by-py*0.35, by+py*0.35], color="#C00000", lw=2.2, zorder=3)

# ---- lncRNA layer: rectangles further out, sponge (dashed) edges to their target miRNA(s) ----
def draw_rect(ax, x, y, label, w, h, facecolor, lw, bold=False, fontsize=8):
    rect = mpatches.FancyBboxPatch((x-w/2, y-h/2), w, h, boxstyle="round,pad=0.05,rounding_size=0.08",
                                     facecolor=facecolor, edgecolor="black", linewidth=lw, zorder=4)
    ax.add_patch(rect)
    ax.text(x, y, label, ha="center", va="center", fontsize=fontsize, fontweight="bold" if bold else "normal", zorder=5)

lncrnas = {
    "AL356488.2": {"targets": ["hsa-miR-654-5p"], "flagship": True, "pos": (9.8, 2.5)},
    "AC020922.3": {"targets": ["hsa-miR-654-5p", "hsa-miR-541-3p"], "flagship": False, "pos": (10.2, -0.5)},
    "AL355987.4": {"targets": ["hsa-miR-541-3p", "hsa-miR-23a-5p"], "flagship": False, "pos": (9.8, -3.2)},
    "LMNTD2-AS1": {"targets": ["hsa-miR-541-3p", "hsa-miR-92a-1-5p"], "flagship": False, "pos": (7.5, -5.8)},
    "NORAD": {"targets": ["hsa-miR-451a"], "flagship": False, "pos": (5.0, -7.3)},
}
for lnc, info in lncrnas.items():
    lx, ly = info["pos"]
    facecolor = "#42A5F5" if info["flagship"] else "#90CAF9"
    lw = 2.8 if info["flagship"] else 1.2
    draw_rect(ax, lx, ly, lnc, 2.0 if len(lnc) > 8 else 1.6, 0.7, facecolor, lw, bold=info["flagship"])
    for target in info["targets"]:
        if target in mi_pos:
            tx, ty = mi_pos[target]
            ax.plot([lx, tx], [ly, ty], color="#1565C0", linewidth=1.6, linestyle="--", zorder=1)

# ---- Legend ----
legend_elements = [
    mpatches.Patch(facecolor="#D32F2F", edgecolor="black", label="MIF"),
    mpatches.Patch(facecolor=tier_colors["full3"], edgecolor="black", label="TF: full 3-species conservation"),
    mpatches.Patch(facecolor=tier_colors["weak3"], edgecolor="black", label="TF: present, sub-threshold"),
    mpatches.Patch(facecolor=tier_colors["heterodimer"], edgecolor="black", label="TF: heterodimer motif only"),
    mpatches.Patch(facecolor=tier_colors["single"], edgecolor="black", label="TF: single-species match"),
    mpatches.Patch(facecolor=tier_colors["none"], edgecolor="black", label="TF: below threshold, all species"),
    mpatches.Patch(facecolor="#F5C242", edgecolor="black", label="miRNA: independently rediscovered"),
    mpatches.Patch(facecolor="#FFE9A8", edgecolor="black", label="miRNA: original panel"),
    mpatches.Patch(facecolor="#42A5F5", edgecolor="black", label="lncRNA: flagship (cross-pipeline convergent)"),
    mpatches.Patch(facecolor="#90CAF9", edgecolor="black", label="lncRNA: other candidate"),
]
ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(-0.32, 1.0), fontsize=8.5, frameon=True,
          title="Legend", title_fontsize=9.5)

ax.annotate("", xy=(-3, 7.8), xytext=(-5, 7.8), arrowprops=dict(arrowstyle="-|>", color="#2E7D32", lw=1.4))
ax.text(-2.8, 7.8, "Activation", va="center", fontsize=8.5, color="#2E7D32")
ax.plot([-3, -5], [7.2, 7.2], color="#C00000", lw=1.4)
ax.plot([-4.9, -4.9], [7.0, 7.4], color="#C00000", lw=2.2)
ax.text(-2.8, 7.2, "Repression", va="center", fontsize=8.5, color="#C00000")
ax.plot([-3, -5], [6.6, 6.6], color="#1565C0", lw=1.6, linestyle="--")
ax.text(-2.8, 6.6, "miRNA sponging", va="center", fontsize=8.5, color="#1565C0")

ax.set_title("MIF Multi-Layer Regulatory Network\nTranscription factors (activation) | MicroRNAs (repression) | lncRNAs (miRNA sponges)",
              fontsize=14, pad=15)

plt.tight_layout()
plt.savefig("MIF_full_multilayer_network.png", dpi=300, bbox_inches="tight")
plt.savefig("MIF_full_multilayer_network.eps", format="eps", bbox_inches="tight")
print("Saved MIF_full_multilayer_network.png/.eps")

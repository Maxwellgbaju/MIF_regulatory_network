import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

fig, ax = plt.subplots(figsize=(13, 11))
ax.set_xlim(-6.5, 6.5)
ax.set_ylim(-6.5, 6.5)
ax.axis("off")

# Center node: MIF -- now red
ax.scatter([0], [0], s=3200, color="#D32F2F", edgecolor="black", linewidth=2, zorder=5)
ax.text(0, 0, "MIF", ha="center", va="center", fontsize=15, fontweight="bold", color="white", zorder=6)

mirnas = ["hsa-miR-320a", "hsa-miR-744-5p", "hsa-miR-146a-5p", "hsa-miR-1246",
          "hsa-miR-769-5p", "hsa-miR-654-5p", "hsa-miR-1228-5p", "hsa-miR-451a"]
no_coverage = {"hsa-miR-1246", "hsa-miR-1228-5p"}
n = len(mirnas)
mirna_pos = {}
for i, m in enumerate(mirnas):
    angle = 2 * np.pi * i / n + np.pi/2  # start from top
    x, y = 2.6 * np.cos(angle), 2.6 * np.sin(angle)
    mirna_pos[m] = (x, y)
    color = "#F5C242" if m not in no_coverage else "#D9D9D9"
    ax.scatter([x], [y], s=1500, color=color, edgecolor="black", linewidth=1.3, zorder=4)
    ax.text(x, y, m, ha="center", va="center", fontsize=8.5, zorder=5)
    ax.plot([0, x], [0, y], color="#5B9BD5", linewidth=1.5, zorder=1)

lncrna_hits = {
    "hsa-miR-320a": ("MIR3142HG", "hub"),
    "hsa-miR-744-5p": ("AL356488.2", "specific"),
    "hsa-miR-146a-5p": ("MALAT1", "hub"),
    "hsa-miR-769-5p": ("MALAT1", "hub"),
    "hsa-miR-654-5p": ("AL356488.2", "specific"),
    "hsa-miR-451a": ("NORAD", "specific"),
}

def draw_rect_node(ax, x, y, label, w, h, facecolor, edgecolor, lw, fontsize, bold=False):
    rect = mpatches.FancyBboxPatch((x - w/2, y - h/2), w, h,
                                     boxstyle="round,pad=0.05,rounding_size=0.08",
                                     facecolor=facecolor, edgecolor=edgecolor, linewidth=lw, zorder=3)
    ax.add_patch(rect)
    ax.text(x, y, label, ha="center", va="center", fontsize=fontsize,
            fontweight="bold" if bold else "normal", zorder=5)

lncrna_unique_pos = {}
placed = set()
for mirna, (lncrna, flag) in lncrna_hits.items():
    mx, my = mirna_pos[mirna]
    if lncrna not in placed:
        angle = np.arctan2(my, mx)
        lx, ly = 5.3 * np.cos(angle), 5.3 * np.sin(angle)
        lncrna_unique_pos[lncrna] = (lx, ly)
        placed.add(lncrna)
    lx, ly = lncrna_unique_pos[lncrna]
    edge_color = "#1565C0" if flag == "specific" else "#AAAAAA"
    ax.plot([mx, lx], [my, ly], color=edge_color, linewidth=2.2 if flag == "specific" else 1.2,
            linestyle="-" if flag == "specific" else "--", zorder=2)

for lncrna, (lx, ly) in lncrna_unique_pos.items():
    is_specific = any(v[0] == lncrna and v[1] == "specific" for v in lncrna_hits.values())
    facecolor = "#5B9BD5" if is_specific else "#EDEDED"
    is_headline = lncrna == "AL356488.2"
    w = 1.9 if len(lncrna) > 6 else 1.5
    draw_rect_node(ax, lx, ly, lncrna, w, 0.75, facecolor, "black",
                   2.5 if is_headline else 1.2, 9, bold=is_headline)

legend_elements = [
    mpatches.Patch(facecolor="#D32F2F", edgecolor="black", label="MIF (the gene of interest)"),
    mpatches.Patch(facecolor="#F5C242", edgecolor="black", label="miRNA -- a linked lncRNA was found"),
    mpatches.Patch(facecolor="#D9D9D9", edgecolor="black", label="miRNA -- no lncRNA data available for this one"),
    mpatches.Patch(facecolor="#5B9BD5", edgecolor="black", label="lncRNA -- specific partner (strong candidate)"),
    mpatches.Patch(facecolor="#EDEDED", edgecolor="black", label="lncRNA -- very common/abundant (weaker evidence)"),
]
ax.legend(handles=legend_elements, loc="lower center", bbox_to_anchor=(0.5, -0.12),
          fontsize=10, frameon=True, ncol=1, title="Legend", title_fontsize=11)

ax.set_title("MIF Regulatory Network: microRNAs and Their Linked lncRNAs",
              fontsize=15, fontweight="bold", pad=20)
ax.text(0, 6.0, "AL356488.2 stands out: two different microRNAs both link to it,\nmaking it the strongest lncRNA candidate in this network.",
        ha="center", fontsize=10.5, style="italic", color="#333333")

plt.tight_layout()
plt.savefig("student_MIF_network_v2.png", dpi=300, bbox_inches="tight")
plt.savefig("student_MIF_network_v2.eps", format="eps", bbox_inches="tight")
print("Saved student_MIF_network_v2.png/.eps")

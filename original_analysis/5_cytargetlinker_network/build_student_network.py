import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(figsize=(12, 10))
ax.set_xlim(-6, 6)
ax.set_ylim(-6, 6)
ax.axis("off")

# Center node: MIF
ax.scatter([0], [0], s=3000, color="#BDBDBD", edgecolor="black", zorder=5)
ax.text(0, 0, "MIF", ha="center", va="center", fontsize=13, fontweight="bold", zorder=6)

# miRNA layer (inner ring)
mirnas = ["hsa-miR-320a", "hsa-miR-744-5p", "hsa-miR-146a-5p", "hsa-miR-1246",
          "hsa-miR-769-5p", "hsa-miR-654-5p", "hsa-miR-1228-5p", "hsa-miR-451a"]
no_coverage = {"hsa-miR-1246", "hsa-miR-1228-5p"}
n = len(mirnas)
mirna_pos = {}
for i, m in enumerate(mirnas):
    angle = 2 * np.pi * i / n
    x, y = 2.5 * np.cos(angle), 2.5 * np.sin(angle)
    mirna_pos[m] = (x, y)
    color = "#F0D264" if m not in no_coverage else "#E0E0E0"
    ax.scatter([x], [y], s=1400, color=color, edgecolor="black", zorder=4)
    ax.text(x, y, m, ha="center", va="center", fontsize=8, zorder=5)
    ax.plot([0, x], [0, y], color="#4CAF50", linewidth=1.3, zorder=1)

# lncRNA layer (outer ring) -- best hit per miRNA, with flags
lncrna_hits = {
    "hsa-miR-320a": ("MIR3142HG", "abundance-caution"),
    "hsa-miR-744-5p": ("AL356488.2", "specific"),
    "hsa-miR-146a-5p": ("MALAT1", "abundance-caution"),
    "hsa-miR-769-5p": ("MALAT1", "abundance-caution"),
    "hsa-miR-654-5p": ("AL356488.2", "specific"),
    "hsa-miR-451a": ("NORAD", "specific"),
}
# Note: hsa-miR-1246 and hsa-miR-1228-5p excluded (no ENCORI lncRNA coverage)

lncrna_unique_pos = {}
placed = set()
for mirna, (lncrna, flag) in lncrna_hits.items():
    mx, my = mirna_pos[mirna]
    if lncrna not in placed:
        # place lncRNA node further out, roughly in direction of its first miRNA
        angle = np.arctan2(my, mx)
        lx, ly = 5.0 * np.cos(angle), 5.0 * np.sin(angle)
        lncrna_unique_pos[lncrna] = (lx, ly)
        placed.add(lncrna)
    lx, ly = lncrna_unique_pos[lncrna]
    edge_color = "#1565C0" if flag == "specific" else "#BDBDBD"
    ax.plot([mx, lx], [my, ly], color=edge_color, linewidth=1.8 if flag == "specific" else 1.0,
            linestyle="-" if flag == "specific" else "--", zorder=2)

for lncrna, (lx, ly) in lncrna_unique_pos.items():
    is_specific = any(v[0] == lncrna and v[1] == "specific" for v in lncrna_hits.values())
    color = "#42A5F5" if is_specific else "#E8E8E8"
    size = 1600 if lncrna == "AL356488.2" else 1200
    ax.scatter([lx], [ly], s=size, color=color, edgecolor="black", zorder=3,
               linewidth=2 if lncrna == "AL356488.2" else 1)
    ax.text(lx, ly, lncrna, ha="center", va="center", fontsize=8,
            fontweight="bold" if lncrna == "AL356488.2" else "normal", zorder=5)

# Legend
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor="#F0D264", edgecolor="black", label="miRNA (ENCORI lncRNA data available)"),
    Patch(facecolor="#E0E0E0", edgecolor="black", label="miRNA (no ENCORI lncRNA coverage)"),
    Patch(facecolor="#42A5F5", edgecolor="black", label="lncRNA: specific, non-ubiquitous target"),
    Patch(facecolor="#E8E8E8", edgecolor="black", label="lncRNA: high-abundance hub (interpret with caution)"),
    Line2D([0], [0], color="#1565C0", linewidth=1.8, label="Specific miRNA-lncRNA link"),
    Line2D([0], [0], color="#BDBDBD", linewidth=1.0, linestyle="--", label="Abundance-driven link (caution)"),
]
ax.legend(handles=legend_elements, loc="upper left", bbox_to_anchor=(-0.05, 1.05), fontsize=8, frameon=True)

ax.set_title("MIF regulatory network: miRNA layer (CyTargetLinker) extended with lncRNA sponges (ENCORI)\nAL356488.2 (bold) is targeted by two independent miRNAs in this network -- the strongest specific candidate",
              fontsize=11)

plt.tight_layout()
plt.savefig("student_MIF_network_v1.png", dpi=300, bbox_inches="tight")
plt.savefig("student_MIF_network_v1.eps", format="eps", bbox_inches="tight")
print("Saved student_MIF_network_v1.png/.eps")

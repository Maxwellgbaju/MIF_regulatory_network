"""Supplementary Figure S2: independent CyTargetLinker microRNA screen and ENCORI lncRNA partners.
Same data as MIF_student_network_v2/build_student_network_v2.py; layout corrected so that every
microRNA-lncRNA link is visible (miR-744-5p and miR-654-5p both link to AL356488.2).
Output: supp_figures/S2_cytargetlinker_network.png and .eps"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "pdf.fonttype": 42, "ps.fonttype": 42})
# microRNAs ordered so that those sharing an lncRNA partner are adjacent
MIRNAS = ["hsa-miR-744-5p", "hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p",
          "hsa-miR-1246", "hsa-miR-769-5p", "hsa-miR-146a-5p", "hsa-miR-320a"]
NO_COVERAGE = {"hsa-miR-1246", "hsa-miR-1228-5p"}          # no lncRNA data in ENCORI
LINKS = [("hsa-miR-744-5p", "AL356488.2", "specific"), ("hsa-miR-654-5p", "AL356488.2", "specific"),
         ("hsa-miR-451a", "NORAD", "specific"), ("hsa-miR-769-5p", "MALAT1", "broad"),
         ("hsa-miR-146a-5p", "MALAT1", "broad"), ("hsa-miR-320a", "MIR3142HG", "broad")]
R_MIR, R_LNC = 2.6, 5.0

fig, ax = plt.subplots(figsize=(8, 8)); ax.set_aspect("equal"); ax.axis("off")
ang = {m: np.pi / 2 - 2 * np.pi * i / len(MIRNAS) + np.pi / len(MIRNAS) for i, m in enumerate(MIRNAS)}
pos = {m: (R_MIR * np.cos(a), R_MIR * np.sin(a)) for m, a in ang.items()}

for m, (x, y) in pos.items():                                  # MIF -> microRNA spokes
    ax.plot([0, x], [0, y], color="#9DBBE0", lw=1.4, zorder=1)
ax.scatter([0], [0], s=3200, color="#D32F2F", edgecolor="black", lw=1.5, zorder=5)
ax.text(0, 0, "MIF", ha="center", va="center", fontsize=13, fontweight="bold", color="white", zorder=6)

lnc_pos, lnc_kind = {}, {}
for lnc in dict.fromkeys(l for _, l, _ in LINKS):              # place each lncRNA between its microRNAs
    a = [ang[m] for m, l, _ in LINKS if l == lnc]
    mean = np.arctan2(np.mean(np.sin(a)), np.mean(np.cos(a)))
    lnc_pos[lnc] = (R_LNC * np.cos(mean), R_LNC * np.sin(mean))
    lnc_kind[lnc] = [k for _, l, k in LINKS if l == lnc][0]

for m, lnc, kind in LINKS:                                     # microRNA -> lncRNA links (curved)
    ax.add_patch(mpatches.FancyArrowPatch(pos[m], lnc_pos[lnc], connectionstyle="arc3,rad=0.12",
                 arrowstyle="-", lw=2.2 if kind == "specific" else 1.2,
                 linestyle="-" if kind == "specific" else "--",
                 color="#1565C0" if kind == "specific" else "#888888", zorder=2))

for m, (x, y) in pos.items():
    ax.scatter([x], [y], s=1700, color="#D9D9D9" if m in NO_COVERAGE else "#F5C242",
               edgecolor="black", lw=1.0, zorder=4)
    ax.text(x, y, m.replace("hsa-", ""), ha="center", va="center", fontsize=8, zorder=5)

for lnc, (x, y) in lnc_pos.items():
    w = 0.32 + 0.16 * len(lnc)
    ax.add_patch(mpatches.FancyBboxPatch((x - w / 2, y - 0.35), w, 0.7, boxstyle="round,pad=0.05,rounding_size=0.08",
                 facecolor="#5B9BD5" if lnc_kind[lnc] == "specific" else "#EDEDED",
                 edgecolor="black", lw=2.0 if lnc == "AL356488.2" else 1.0, zorder=3))
    ax.text(x, y, lnc, ha="center", va="center", fontsize=9,
            fontweight="bold" if lnc == "AL356488.2" else "normal", zorder=5)

lim = R_LNC + 1.4; ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
handles = [mpatches.Patch(facecolor="#F5C242", edgecolor="black", label="microRNA with ENCORI lncRNA partner"),
           mpatches.Patch(facecolor="#D9D9D9", edgecolor="black", label="microRNA without ENCORI lncRNA data"),
           mpatches.Patch(facecolor="#5B9BD5", edgecolor="black", label="lncRNA, specific partner"),
           mpatches.Patch(facecolor="#EDEDED", edgecolor="black", label="lncRNA, broadly expressed"),
           Line2D([0], [0], color="#9DBBE0", lw=1.4, label="Predicted microRNA targeting of MIF (CyTargetLinker)"),
           Line2D([0], [0], color="#1565C0", lw=2.2, label="CLIP-supported interaction, specific lncRNA (ENCORI)"),
           Line2D([0], [0], color="#888888", lw=1.2, ls="--", label="CLIP-supported interaction, broadly expressed lncRNA")]
ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=2, fontsize=7.5, frameon=False)
os.makedirs("supp_figures", exist_ok=True)
fig.savefig("supp_figures/S2_cytargetlinker_network.png", dpi=600, bbox_inches="tight")
fig.savefig("supp_figures/S2_cytargetlinker_network.eps", format="eps", bbox_inches="tight")
print("Saved supp_figures/S2_cytargetlinker_network.png and .eps")

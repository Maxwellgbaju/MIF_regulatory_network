"""Figure 8: microRNA-lncRNA predictions and their expression test. Input: results/rev4_mirna_cerna_summary.csv (rev4).
miR-608 is excluded: rev4 matched it to the wrong microRNA (hsa-miR-6087)."""
import numpy as np, pandas as pd
from figstyle import plt, save, tag, need, RED, BLUE, GREY
from matplotlib.patches import FancyBboxPatch

s = pd.read_csv(need("results/rev4_mirna_cerna_summary.csv"))
s = s[~s.a.isin(["hsa-miR-608"])]
fig = plt.figure(figsize=(7.4, 6.0))
gs = fig.add_gridspec(2, 2, height_ratios=[0.7, 1.3], hspace=0.5, wspace=0.75)

a = fig.add_subplot(gs[0, :]); a.axis("off"); a.set_xlim(0, 100); a.set_ylim(0, 100)
def box(x, y, w, h, title, body, fc):
    a.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.5,rounding_size=1.5", fc=fc, ec="#333", lw=0.8))
    a.text(x + w / 2, y + h - 5, title, ha="center", va="top", fontsize=7.2, fontweight="bold")
    a.text(x + w / 2, y + h - 17, body, ha="center", va="top", fontsize=6.2, linespacing=1.3)
box(1, 8, 27, 86, "Screen 1 (multiMiR)", "miR-92a-1-5p, miR-608,\nmiR-541-3p, miR-363-5p,\nmiR-23a-5p, miR-23b-5p,\nmiR-1537-5p", "#E8F0FA")
box(36.5, 8, 27, 86, "Shared by both screens", "miR-654-5p\nmiR-451a\nmiR-1228-5p\n\nShared lncRNA partner: AL356488.2\n(linked to miR-654-5p in both screens)", "#FFF4D6")
box(72, 8, 27, 86, "Screen 2 (CyTargetLinker)", "miR-744-5p, miR-320a,\nmiR-146a-5p, miR-1246,\nmiR-769-5p", "#E9F5EC")
for x0, x1 in [(28.6, 36), (71.4, 64)]:
    a.annotate("", xy=(x1, 50), xytext=(x0, 50), arrowprops=dict(arrowstyle="-|>", lw=1.0, color="#333"))
a.text(-0.03, 1.02, "a", transform=a.transAxes, fontsize=11, fontweight="bold")

def panel(ax, test, title, expected, letter):
    d = s[s.test == test].sort_values("median_rho").reset_index(drop=True); y = np.arange(len(d))
    ax.barh(y, d.sig_expected_dir, color=RED if expected == "positive" else BLUE, height=0.6, label=f"expected direction ({expected})")
    ax.barh(y, -d.sig_opposite_dir, color=GREY, height=0.6, label="opposite direction")
    labels = d.a.str.replace("hsa-", "", regex=False) if test != "miRNA_vs_lncRNA" else (d.a.str.replace("hsa-", "", regex=False) + " / " + d.b)
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=6.2); ax.axvline(0, color="k", lw=0.6)
    lim = max(5, int(max(d.sig_expected_dir.max(), d.sig_opposite_dir.max())) + 3)
    ax.set_xlim(-lim, lim + 6)
    for yi, m in zip(y, d.median_rho): ax.text(lim + 5.5, yi, f"{m:+.2f}", va="center", ha="right", fontsize=5.8, color="#444")
    ax.set_xlabel("Cancer types with FDR < 0.05\n(opposite | expected direction)"); ax.legend(fontsize=5.6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.3))
    ax.set_title(title, fontsize=8); tag(ax, letter)
panel(fig.add_subplot(gs[1, 0]), "miRNA_vs_MIF", "microRNA vs MIF\n(sponge model expects negative)", "negative", "b")
panel(fig.add_subplot(gs[1, 1]), "lncRNA_vs_MIF", "lncRNA vs MIF\n(sponge model expects positive)", "positive", "c")
save(fig, "Figure8_mirna_lncrna")

"""Figure 3: co-expression of candidates with MIF. Input: results/rev1_MIF_coexpression_all.csv (rev1)."""
import numpy as np, pandas as pd
from figstyle import plt, save, tag, need, RED, BLUE, GREY

r = pd.read_csv(need("results/rev1_MIF_coexpression_all.csv"))
s = (r.assign(pos=(r.rho > 0) & (r.FDR < 0.05), neg=(r.rho < 0) & (r.FDR < 0.05))
       .groupby("gene").agg(pos=("pos", "sum"), neg=("neg", "sum"), med=("rho", "median")).sort_values("med"))
t = r[r.gene == "TFAP2A"].sort_values("rho", ascending=False)

fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 4.8), gridspec_kw=dict(width_ratios=[1.05, 1], wspace=0.6))
y = np.arange(len(s))
a.barh(y, s.pos, color=RED, height=0.65, label="significantly positive")
a.barh(y, -s.neg, color=BLUE, height=0.65, label="significantly negative")
a.set_yticks(y); a.set_yticklabels(s.index, fontsize=7); a.axvline(0, color="k", lw=0.6)
lim = 5 * int(np.ceil(max(s.pos.max(), s.neg.max()) / 5)) + 5
a.set_xlim(-lim, lim); ticks = np.arange(-lim + 5, lim, 10); a.set_xticks(ticks); a.set_xticklabels([str(abs(int(x))) for x in ticks])
a.set_xlabel("Cancer types with FDR < 0.05\n(negative | positive association with MIF)")
for yi, m in zip(y, s.med): a.text(lim - 0.5, yi, f"{m:+.2f}", va="center", ha="right", fontsize=6.2, color="#444")
a.text(lim - 0.5, len(s) - 0.3, "median rho", ha="right", fontsize=6.2, color="#444")
a.legend(fontsize=6.2, frameon=False, loc="upper left")
a.set_title(f"Candidates vs MIF, {r.cancer.nunique()} cancer types", fontsize=8.3); tag(a, "a")
yb = np.arange(len(t))
cols = [RED if (f < .05 and x > 0) else BLUE if (f < .05 and x < 0) else GREY for x, f in zip(t.rho, t.FDR)]
b.hlines(yb, 0, t.rho, color=cols, lw=1.4); b.scatter(t.rho, yb, color=cols, s=14, zorder=3)
b.set_yticks(yb); b.set_yticklabels(t.cancer, fontsize=5.6); b.invert_yaxis(); b.axvline(0, color="k", lw=0.6)
b.set_xlabel("Spearman rho, TFAP2A vs MIF"); b.set_title("TFAP2A: context-dependent association", fontsize=8.3); tag(b, "b")
save(fig, "Figure3_coexpression")

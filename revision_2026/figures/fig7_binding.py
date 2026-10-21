"""Figure 7: ChIP-Atlas binding at the MIF promoter. Input: results/rev5_chipatlas_MIF_promoter.csv (rev5)."""
import numpy as np, pandas as pd
from figstyle import plt, save, tag, need, RED

ca = pd.read_csv(need("results/rev5_chipatlas_MIF_promoter.csv"))
ca = ca[ca.experiments_total > 0].copy()
ca["pct_bound"] = 100 * ca.experiments_binding_MIF.fillna(0) / ca.experiments_total
ca["perc"] = ca.MIF_score_percentile_among_targets.fillna(-1)
ca = ca.sort_values("perc", ascending=False).reset_index(drop=True)
hi = {"HIF1A", "TFAP2A"}; cols = [RED if t in hi else "#4C72B0" for t in ca.TF]; y = np.arange(len(ca))
fig, (a, b) = plt.subplots(1, 2, figsize=(7.2, 4.4), sharey=True, gridspec_kw=dict(wspace=0.08))
a.barh(y, ca.pct_bound, color=cols, height=0.65)
for yi, r in ca.iterrows(): a.text(r.pct_bound + 1.5, yi, f"{int(r.experiments_binding_MIF if pd.notna(r.experiments_binding_MIF) else 0)}/{int(r.experiments_total)}", va="center", fontsize=6)
a.set_yticks(y); a.set_yticklabels(ca.TF, fontsize=7); a.invert_yaxis(); a.set_xlim(0, max(85, ca.pct_bound.max() + 15))
a.set_xlabel("% of ChIP-seq experiments with a peak\nwithin 1 kb of the MIF TSS"); tag(a, "a")
b.barh(y, ca.perc.clip(lower=0), color=cols, height=0.65); b.axvline(90, color="k", lw=0.6, ls="--")
for yi, r in ca.iterrows():
    if r.perc < 0: b.text(2, yi, "no peak", va="center", fontsize=6)
b.set_xlim(0, 105); b.set_xlabel("MIF binding score: percentile\namong all genes bound by the TF"); tag(b, "b")
fig.suptitle("Direct binding at the MIF promoter (ChIP-Atlas, hg38)", fontsize=9, y=0.99)
save(fig, "Figure7_binding")

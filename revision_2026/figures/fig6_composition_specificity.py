"""Figure 6: composition and specificity. Inputs: results/rev3_tf_activity_vs_MIF.csv, results/rev6_hypoxia_vs_MIF.csv,
results/rev7_hallmark_specificity_summary.csv, results/rev8_program_specificity.csv."""
import numpy as np, pandas as pd
from figstyle import plt, save, tag, need, RED, GREY

r3 = pd.read_csv(need("results/rev3_tf_activity_vs_MIF.csv"))
r6 = pd.read_csv(need("results/rev6_hypoxia_vs_MIF.csv"))
r7 = pd.read_csv(need("results/rev7_hallmark_specificity_summary.csv"), index_col=0)
r8 = pd.read_csv(need("results/rev8_program_specificity.csv"))
fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.4), gridspec_kw=dict(hspace=0.62, wspace=0.75))

a = axs[0, 0]
s3 = r3.groupby("TF").agg(m=("rho_mRNA_vs_MIF", "median"), act=("rho_activity_vs_MIF", "median")).sort_values("act", ascending=False)
nsig = int((r3.FDR_activity_perm < 0.05).sum()); y = np.arange(len(s3))
for yi, (_, r) in zip(y, s3.iterrows()): a.plot([r.m, r.act], [yi, yi], color="#CCCCCC", lw=0.8, zorder=0)
a.scatter(s3.m, y, color="#7F7F7F", s=16, label="TF mRNA"); a.scatter(s3.act, y, color="#8E44AD", s=16, marker="D", label="TF activity (targets)")
a.set_yticks(y); a.set_yticklabels(s3.index, fontsize=6.6); a.axvline(0, color="k", lw=0.6)
a.set_xlabel(f"Median rho with MIF ({r3.cancer.nunique()} cancer types)"); a.legend(fontsize=6, frameon=False, loc="upper right")
a.set_title(f"TF activity vs random target sets\n({nsig} significant TF-cancer pairs)", fontsize=8); tag(a, "a")

b = axs[0, 1]; d6 = r6.sort_values("partial_rho_adj_stroma_immune", ascending=False).reset_index(drop=True); y = np.arange(len(d6))
hit = (d6.FDR < 0.05) & (d6.partial_rho_adj_stroma_immune > d6.null_partial_median)
for yi, r in d6.iterrows(): b.plot([r.null_partial_median, r.partial_rho_adj_stroma_immune], [yi, yi], color="#DDDDDD", lw=0.7, zorder=0)
b.scatter(d6.null_partial_median, y, color=GREY, s=10, label="random gene sets (median)")
b.scatter(d6.partial_rho_adj_stroma_immune, y, color=[RED if h else "#555555" for h in hit], s=12, label="Hallmark hypoxia")
b.axvline(0, color="k", lw=0.6); b.set_yticks([]); b.invert_yaxis(); b.set_ylabel(f"{len(d6)} cancer types")
b.set_xlabel("Composition-adjusted rho with MIF"); b.legend(fontsize=6, frameon=False, loc="upper left")
b.set_title(f"Hypoxia vs random sets\n(red: exceeds random, FDR < 0.05; {int(hit.sum())} of {len(d6)})", fontsize=8); tag(b, "b")

c = axs[1, 0]; top = r7.sort_values("median_rank").head(10); y = np.arange(len(top))
c.barh(y, top.median_rank, color=[RED if n == "Hypoxia" else "#4C72B0" for n in top.index], height=0.65)
for yi, (n, r) in zip(y, top.iterrows()): c.text(r.median_rank + 0.3, yi, f"top-5 in {int(r.cancers_in_top5)}", va="center", fontsize=6)
pos = list(r7.sort_values("median_rank").index).index("Hypoxia") + 1
c.set_yticks(y); c.set_yticklabels(top.index, fontsize=6.6); c.invert_yaxis(); c.set_xlim(0, top.median_rank.max() + 7)
c.set_xlabel(f"Median rank among {len(r7)} Hallmark programmes"); c.set_title(f"Programmes most associated with MIF\n(hypoxia ranks {pos})", fontsize=8); tag(c, "c")

d = axs[1, 1]
s8 = r8.groupby("program").agg(pct=("MIF_percentile_vs_controls", "median"),
                               above=("MIF_percentile_vs_controls", lambda x: int((x >= 90).sum()))).sort_values("pct")
y = np.arange(len(s8))
d.barh(y, s8.pct, color="#4C72B0", height=0.6); d.axvline(90, color=RED, lw=0.9, ls="--"); d.axvline(50, color="k", lw=0.6, ls=":")
d.set_yticks(y); d.set_yticklabels(s8.index, fontsize=6.6); d.set_xlim(0, 100)
d.set_xlabel("MIF percentile vs 100 expression-\nmatched genes (median over cancers)")
d.set_title(f"MIF vs expression-matched genes\n({int(s8.above.max())} cancers above 90th, any programme)", fontsize=8); tag(d, "d")
save(fig, "Figure6_composition_specificity")

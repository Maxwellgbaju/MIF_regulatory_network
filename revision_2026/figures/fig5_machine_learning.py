"""Figure 5: cross-database machine learning (kidney clear cell carcinoma).
Inputs: results/ml_stability_selection.csv, results/ml_cindex.csv, results/ml_rsf_importance.csv (make_ml_values.py)."""
import numpy as np, pandas as pd
from figstyle import plt, save, tag, need, RED

st = pd.read_csv(need("results/ml_stability_selection.csv"))
ci = pd.read_csv(need("results/ml_cindex.csv"))
rf = pd.read_csv(need("results/ml_rsf_importance.csv"))
fig, axs = plt.subplots(1, 3, figsize=(7.6, 4.2), gridspec_kw=dict(wspace=0.85, width_ratios=[1.1, 1.1, 1]))

a = axs[0]
st = st.assign(order=st[["stability_cbioportal", "stability_xena"]].max(axis=1)).sort_values("order")
y = np.arange(len(st))
for yi, (_, r) in zip(y, st.iterrows()):
    if pd.notna(r.stability_cbioportal) and pd.notna(r.stability_xena):
        a.plot([r.stability_cbioportal, r.stability_xena], [yi, yi], color="#DDDDDD", lw=1, zorder=0)
a.scatter(st.stability_cbioportal, y, color="#7F7F7F", s=16, label="cBioPortal (original quantification)")
a.scatter(st.stability_xena, y, color="#009E73", s=16, marker="D", label="Xena (Toil)")
a.axvline(80, color="k", lw=0.6, ls="--")
a.set_yticks(y); a.set_yticklabels(st.gene, fontsize=6.6); a.set_xlim(-3, 103)
a.set_xlabel("LASSO stability selection (%)"); a.legend(fontsize=5.6, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22))
a.set_title("Selected features differ\nbetween databases", fontsize=8); tag(a, "a")
for yi, g in zip(y, st.gene):
    if g == "RELA": a.get_yticklabels()[list(st.gene).index(g)].set_color(RED)

b = axs[1]
models = ["Multi-gene LASSO", "RELA alone", "MIF alone", "Random survival forest"]
for k, (ds, col) in enumerate(zip(["cBioPortal (n=510)", "Xena (n=530)"], ["#7F7F7F", "#009E73"])):
    sub = ci[ci.dataset == ds].set_index("model")
    xs = [i + (k - 0.5) * 0.36 for i, m in enumerate(models) if m in sub.index]
    vals = [sub.loc[m, "c_index"] for m in models if m in sub.index]
    sds = [sub.loc[m, "sd"] for m in models if m in sub.index]
    b.bar(xs, vals, 0.34, yerr=sds, color=col, label=ds.split(" (")[0], error_kw=dict(lw=0.7, capsize=2))
b.axhline(0.5, color="k", lw=0.6, ls=":"); b.set_ylim(0.45, 0.72)
b.set_xticks(range(len(models))); b.set_xticklabels(["LASSO\n(multi-\ngene)", "RELA\nonly", "MIF\nonly", "RSF"], fontsize=5.8)
b.set_ylabel("Cross-validated C-index (mean, SD)"); b.legend(fontsize=6, frameon=False, loc="upper left")
b.set_title("Discrimination by model", fontsize=8); tag(b, "b")

c = axs[2]; rf = rf.sort_values("importance")
c.barh(rf.gene, rf.importance, xerr=rf.sd, color=[RED if g == "TFAP2A" else "#4C72B0" for g in rf.gene], height=0.65,
       error_kw=dict(lw=0.6, capsize=1.5))
c.tick_params(axis="y", labelsize=6.6); c.set_xlabel("Permutation importance\n(mean C-index drop)")
c.set_title("Random survival forest\n(cBioPortal)", fontsize=8); tag(c, "c")
save(fig, "Figure5_machine_learning")

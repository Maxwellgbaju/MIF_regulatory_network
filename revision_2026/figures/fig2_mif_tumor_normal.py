"""Figure 2: MIF tumour vs normal. Recomputes the comparison from the Toil data (same grouping as the
September analysis: TCGA primary tumours vs TCGA solid tissue normal of the same disease plus GTEx normal
of the same primary site; >=10 samples per group), saves results/fig2_MIF_tumor_vs_normal.csv
(Supplementary Table S1), then plots it."""
import numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
from mifrev_utils import load_meta, resolve, fetch, bh
from figstyle import plt, save

meta = load_meta()
tum, rows = {}, []
for s, m in meta.items():
    if m.get("_study") == "TCGA" and m.get("_sample_type") == "Primary Tumor":
        tum.setdefault(m.get("primary disease or tissue"), []).append(s)
groups = {}
for dis, t in tum.items():
    if not dis or len(t) < 10:
        continue
    sites = {meta[s].get("_primary_site") for s in t if meta[s].get("_primary_site")}
    if not sites:
        continue
    site = sorted(sites)[0]
    n_tcga = [s for s, m in meta.items() if m.get("_study") == "TCGA" and m.get("_sample_type") == "Solid Tissue Normal"
              and m.get("primary disease or tissue") == dis]
    n_gtex = [s for s, m in meta.items() if m.get("_study") == "GTEX" and m.get("_sample_type") == "Normal Tissue"
              and m.get("_primary_site") == site]
    if len(n_tcga) + len(n_gtex) >= 10:
        groups[dis] = (t, n_tcga, n_gtex, site)
samples = sorted({s for t, a, b, _ in groups.values() for s in t + a + b})
ids, _ = resolve(["MIF"])
mif = fetch(ids, samples, "fig2_tvn")["MIF"]
for dis, (t, a, b, site) in groups.items():
    tv, nv = mif.reindex(t).dropna(), mif.reindex(a + b).dropna()
    p = mannwhitneyu(tv, nv, alternative="two-sided").pvalue
    rows.append({"cancer": dis, "primary_site": site, "n_tumor": len(tv), "n_normal_TCGA": len(a), "n_normal_GTEx": len(b),
                 "median_tumor": tv.median(), "median_normal": nv.median(),
                 "log2FC_median": tv.median() - nv.median(), "p": p})
res = pd.DataFrame(rows)
res["FDR"] = bh(res.p)
res = res.sort_values("log2FC_median")
res.to_csv("results/fig2_MIF_tumor_vs_normal.csv", index=False)
up = ((res.FDR < 0.05) & (res.log2FC_median > 0)).sum()
print(f"{len(res)} cancer types tested; MIF significantly up in {up}; down in {((res.FDR < 0.05) & (res.log2FC_median < 0)).sum()}")

fig, ax = plt.subplots(figsize=(5.2, max(4.0, 0.22 * len(res))))
lp = -np.log10(res.FDR.clip(lower=1e-300)).values
cmap = plt.cm.viridis
ax.barh(range(len(res)), res.log2FC_median, color=cmap(np.clip(lp / 110, 0, 1)), height=0.72)
ax.set_yticks(range(len(res))); ax.set_yticklabels(res.cancer, fontsize=6.5)
ax.axvline(0, color="k", lw=0.6); ax.set_xlabel("MIF log$_2$ fold change (median), tumour vs normal")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 110))
fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.02).set_label("-log$_{10}$ FDR (capped at 110)", fontsize=7)
ax.set_title(f"MIF upregulated in {up} of {len(res)} cancer types (FDR < 0.05)", fontsize=8.5)
ax.text(0.98, 0.03, "MIF somatic mutation: 1 of 2,688 tumours (0.04%)\nin six TCGA cohorts (Supplementary Fig. S3)",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=6.5, bbox=dict(fc="white", ec="#999", lw=0.5))
save(fig, "Figure2_MIF_tumor_vs_normal")

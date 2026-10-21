from mifrev_utils import *
from scipy.stats import spearmanr
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

GENES = ["MIF", "SP1", "HIF1A", "RELA", "EGR1", "NR1H2", "RXRA", "IRF1", "KLF4",
         "IRF2", "NFIL3", "ESR1", "TP53", "TFAP2A", "NORAD", "AL356488.2"]
meta = load_meta()
groups = tumor_groups(meta, 30)
samples = sorted({s for v in groups.values() for s in v})
ids, missing = resolve(GENES)
print(f"{len(groups)} cancer types, {len(samples)} tumors; unresolved genes: {missing}")
expr = fetch(ids, samples, "panel")

rows = []
for dis, ss in groups.items():
    sub = expr.loc[ss]
    for g in ids:
        if g == "MIF":
            continue
        d = sub[["MIF", g]].dropna()
        if len(d) < 30 or d[g].nunique() < 10:
            continue
        rho, p = spearmanr(d["MIF"], d[g])
        rows.append({"cancer": dis, "gene": g, "n": len(d), "rho": rho, "p": p})
res = pd.DataFrame(rows)
res["FDR"] = bh(res["p"])
res.to_csv("results/rev1_MIF_coexpression_all.csv", index=False)

summ = (res.assign(pos=(res.rho > 0) & (res.FDR < 0.05), neg=(res.rho < 0) & (res.FDR < 0.05))
           .groupby("gene").agg(cancers_tested=("cancer", "count"), sig_positive=("pos", "sum"),
                                sig_negative=("neg", "sum"), median_rho=("rho", "median"))
           .sort_values("median_rho", ascending=False))
summ.to_csv("results/rev1_MIF_coexpression_summary.csv")
print("\n=== Co-expression with MIF across tumor types (Spearman, BH-FDR over all tests) ===")
print(summ.round(3).to_string())
print("\n=== TFAP2A by cancer type ===")
print(res[res.gene == "TFAP2A"].sort_values("rho")[["cancer", "n", "rho", "FDR"]].round(4).to_string(index=False))

mat = res.pivot(index="gene", columns="cancer", values="rho").loc[summ.index]
fdr = res.pivot(index="gene", columns="cancer", values="FDR").loc[summ.index]
fig, ax = plt.subplots(figsize=(max(10, 0.45 * mat.shape[1]), 6.5))
im = ax.imshow(mat.values, cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        if fdr.values[i, j] < 0.05:
            ax.text(j, i, "*", ha="center", va="center", fontsize=8)
ax.set_xticks(range(mat.shape[1])); ax.set_xticklabels(mat.columns, rotation=60, ha="right", fontsize=7)
ax.set_yticks(range(mat.shape[0])); ax.set_yticklabels(mat.index, fontsize=8)
plt.colorbar(im, ax=ax, label="Spearman rho with MIF", shrink=0.7)
ax.set_title("Co-expression of candidate regulators with MIF in tumors (* FDR < 0.05)", fontsize=10)
plt.tight_layout()
plt.savefig("results/rev1_MIF_coexpression_heatmap.png", dpi=600)
plt.savefig("results/rev1_MIF_coexpression_heatmap.eps")
print("\nSaved results/rev1_*")

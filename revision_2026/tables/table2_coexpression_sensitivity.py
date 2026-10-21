"""Table 2: sensitivity of MIF-candidate co-expression to processing pipeline and analysis parameters.
Eight factors x seven cohorts x two pipelines (Toil re-quantification; original TCGA RSEM) x four settings:
(1) Spearman, all tumours; (2) Pearson on log expression; (3) Spearman adjusted for stromal and immune scores;
(4) Spearman in tumour-rich samples (lower half of the combined stromal+immune score).
Outputs: tables/Table2_long.csv (every test) and tables/Table2_coexpression_sensitivity.csv (summary)."""
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
from tab_common import COHORTS, STROMA, IMMUNE, legacy_expr, toil_expr, composition_score, partial_spearman, bh, show

GENES = ["RELA", "SP1", "HIF1A", "TFAP2A", "NR1H2", "MYC", "EGR1", "KLF4"]
ALL = ["MIF"] + GENES + STROMA + IMMUNE
toil, bars = toil_expr(ALL, "table2")
rows = []
for dis, code in COHORTS.items():
    print(f"{code} ...", flush=True)
    sources = {"Toil": toil.reindex(bars[dis])}
    try:
        sources["Original"] = legacy_expr(code, ALL).reindex(bars[dis]).dropna(how="all")
    except Exception as ex:
        print(f"  original data for {code} failed: {ex}")
    for src, e in sources.items():
        e = e.dropna(subset=["MIF"])
        st, im = composition_score(e)
        comp = (st + im) / 2
        rich = comp <= comp.median()
        for g in GENES:
            if g not in e:
                continue
            d = pd.concat([e["MIF"], e[g], st, im, rich.rename("rich")], axis=1, keys=["MIF", "g", "st", "im", "rich"]).dropna()
            if len(d) < 30:
                continue
            settings = {
                "Spearman, all": spearmanr(d.MIF, d.g),
                "Pearson, all": pearsonr(d.MIF, d.g),
                "Spearman, composition-adjusted": partial_spearman(d.g, d.MIF, [d.st, d.im]),
                "Spearman, tumour-rich half": spearmanr(d[d.rich.astype(bool)].MIF, d[d.rich.astype(bool)].g)}
            for k, (rho, p) in settings.items():
                rows.append({"cohort": code, "gene": g, "pipeline": src, "setting": k, "n": len(d) if "rich" not in k else int(d.rich.sum()),
                             "rho": float(rho), "p": float(p)})
L = pd.DataFrame(rows)
L["FDR"] = bh(L.p)
L.to_csv("tables/Table2_long.csv", index=False)

def summ(x):
    pos = int(((x.rho > 0) & (x.FDR < 0.05)).sum()); neg = int(((x.rho < 0) & (x.FDR < 0.05)).sum())
    return f"{x.rho.median():+.2f} ({pos}+/{neg}-)"
S = L.groupby(["gene", "pipeline", "setting"]).apply(summ).unstack("setting")
S = S[["Spearman, all", "Pearson, all", "Spearman, composition-adjusted", "Spearman, tumour-rich half"]]
S.to_csv("tables/Table2_coexpression_sensitivity.csv")
print("\nCells: median rho across 7 cohorts (number of cohorts significantly positive / negative, FDR < 0.05)")
show("TABLE 2", S)

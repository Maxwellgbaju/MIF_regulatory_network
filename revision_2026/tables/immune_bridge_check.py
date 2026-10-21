"""Immune bridge check: MIF vs PDCD1, FOXP3, CD163, CD8A, CD8B in seven TCGA cohorts, Toil vs original TCGA RSEM.
Settings: Spearman; partial Spearman adjusted for PTPRC; adjusted for PTPRC and stromal score.
Outputs: tables/immune_bridge_long.csv, tables/immune_bridge_summary.csv, tables/immune_bridge_agreement.csv"""
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from tab_common import COHORTS, STROMA, legacy_expr, toil_expr, partial_spearman, zscore, bh, show

MARKERS = ["PDCD1", "FOXP3", "CD163", "CD8A", "CD8B"]
ALL = ["MIF", "PTPRC"] + MARKERS + STROMA
toil, bars = toil_expr(ALL, "immune_bridge")
rows, agree = [], []
for dis, code in COHORTS.items():
    print(f"{code} ...", flush=True)
    sources = {"Toil": toil.reindex(bars[dis])}
    try:
        sources["Original"] = legacy_expr(code, ALL).reindex(bars[dis])
    except Exception as ex:
        print(f"  original data for {code} failed: {ex}")
    if "Original" in sources:
        for g in ["MIF", "PTPRC"] + MARKERS:
            d = pd.concat([sources["Toil"][g], sources["Original"][g]], axis=1).dropna()
            if len(d) >= 30:
                agree.append({"cohort": code, "gene": g, "rho_toil_vs_original": spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0]})
    for src, e in sources.items():
        e = e.dropna(subset=["MIF", "PTPRC"])
        z = zscore(e)
        stroma = z[[g for g in STROMA if g in z]].mean(axis=1)
        for m in MARKERS:
            if m not in e:
                continue
            d = pd.concat([e["MIF"], e[m], e["PTPRC"], stroma], axis=1, keys=["MIF", "m", "PTPRC", "st"]).dropna()
            if len(d) < 30 or d.m.nunique() < 10:
                continue
            for setting, (rho, p) in {"Spearman": spearmanr(d.MIF, d.m),
                                      "Partial | PTPRC": partial_spearman(d.m, d.MIF, [d.PTPRC]),
                                      "Partial | PTPRC + stroma": partial_spearman(d.m, d.MIF, [d.PTPRC, d.st])}.items():
                rows.append({"cohort": code, "marker": m, "pipeline": src, "setting": setting,
                             "n": len(d), "rho": float(rho), "p": float(p)})
L = pd.DataFrame(rows)
L["FDR"] = bh(L.p)
L.to_csv("tables/immune_bridge_long.csv", index=False)
def summ(x):
    pos = int(((x.rho > 0) & (x.FDR < 0.05)).sum()); neg = int(((x.rho < 0) & (x.FDR < 0.05)).sum())
    return f"{x.rho.median():+.2f} ({pos}+/{neg}-)"
S = L.groupby(["marker", "pipeline", "setting"]).apply(summ).unstack("setting")[["Spearman", "Partial | PTPRC", "Partial | PTPRC + stroma"]]
S.to_csv("tables/immune_bridge_summary.csv")
print("Cells: median rho across cohorts (cohorts significantly positive/negative, FDR < 0.05)")
show("IMMUNE SUMMARY", S)
W = L[L.setting == "Partial | PTPRC"].pivot_table(index=["marker", "cohort"], columns="pipeline", values="rho").round(3)
show("IMMUNE PTPRC-ADJUSTED BY COHORT", W)
A = pd.DataFrame(agree)
A.to_csv("tables/immune_bridge_agreement.csv", index=False)
show("IMMUNE PIPELINE AGREEMENT", A.pivot(index="gene", columns="cohort", values="rho_toil_vs_original").round(3))

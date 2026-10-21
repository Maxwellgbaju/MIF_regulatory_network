"""Table 3: sensitivity of survival associations to pipeline, endpoint and model.
Genes: MIF, RELA, SP1, HIF1A, TFAP2A. Seven cohorts. Pipelines: Toil re-quantification and original TCGA RSEM.
Endpoints: overall (OS), disease-specific (DSS) and progression-free interval (PFI), from the TCGA Pan-Cancer Clinical
Data Resource (Xena 'TCGA_survival_data'). Models: Cox per SD of expression; Cox high vs low (median split).
Outputs: tables/Table3_long.csv (every model) and tables/Table3_survival_sensitivity.csv (formatted)."""
import numpy as np, pandas as pd
from lifelines import CoxPHFitter
from mifrev_utils import TOIL, load_meta, tumor_groups
from tab_common import xena, COHORTS, legacy_expr, toil_expr, bh, stars, show

GENES = ["MIF", "RELA", "SP1", "HIF1A", "TFAP2A"]
meta = load_meta(); groups = tumor_groups(meta, 30)
all_t = sorted({s for v in groups.values() for s in v})
fields = xena.dataset_field(TOIL, "TCGA_survival_data")
ENDS = [e for e in ["OS", "DSS", "PFI"] if e in fields and f"{e}.time" in fields]
print("Endpoints available:", ENDS)
sv = xena.dataset_probe_values(TOIL, "TCGA_survival_data", all_t, sum([[e, f"{e}.time"] for e in ENDS], []))
surv = pd.DataFrame({name: pd.to_numeric(pd.Series(v, index=all_t), errors="coerce")
                     for name, v in zip(sum([[e, f"{e}.time"] for e in ENDS], []), sv[1])})
surv.index = [s[:15] for s in surv.index]; surv = surv[~surv.index.duplicated()]

toil, bars = toil_expr(GENES, "table3")
rows = []
for dis, code in COHORTS.items():
    print(f"{code} ...", flush=True)
    sources = {"Toil": toil.reindex(bars[dis])}
    try:
        sources["Original"] = legacy_expr(code, GENES).reindex(bars[dis])
    except Exception as ex:
        print(f"  original data for {code} failed: {ex}")
    for src, e in sources.items():
        for g in GENES:
            for end in ENDS:
                d = pd.concat([e[g].rename("x"), surv.reindex(e.index)[[end, f"{end}.time"]]], axis=1).dropna()
                d.columns = ["x", "event", "time"]; d = d[d.time > 0]
                if len(d) < 50 or d.event.sum() < 10:
                    continue
                for model in ["per SD", "high vs low"]:
                    dd = d.copy()
                    dd["z"] = (dd.x - dd.x.mean()) / dd.x.std() if model == "per SD" else (dd.x > dd.x.median()).astype(int)
                    try:
                        s = CoxPHFitter().fit(dd[["z", "time", "event"]], "time", "event").summary.loc["z"]
                    except Exception as ex:
                        print(f"  {code} {g} {end} {model}: {ex}"); continue
                    rows.append({"cohort": code, "gene": g, "pipeline": src, "endpoint": end, "model": model,
                                 "n": len(dd), "events": int(dd.event.sum()), "HR": s["exp(coef)"],
                                 "CI_low": s["exp(coef) lower 95%"], "CI_high": s["exp(coef) upper 95%"], "p": s["p"]})
L = pd.DataFrame(rows)
L["FDR"] = bh(L.p)
L.to_csv("tables/Table3_long.csv", index=False)

P = L[L.model == "per SD"].copy()
P["cell"] = [f"{h:.2f}{stars(f)}" for h, f in zip(P.HR, P.FDR)]
W = P.pivot_table(index=["gene", "cohort"], columns=["endpoint", "pipeline"], values="cell", aggfunc="first")
W = W.reindex(columns=[(e, p) for e in ENDS for p in ["Toil", "Original"] if (e, p) in W.columns])
W.columns = [f"{e} {p}" for e, p in W.columns]
W.to_csv("tables/Table3_survival_sensitivity.csv")
print("\nHR per SD; * FDR<0.05, ** <0.01, *** <0.001 (corrected across all models)")
show("TABLE 3", W)
agree = L.pivot_table(index=["cohort", "gene", "endpoint", "model"], columns="pipeline", values="FDR")
if {"Toil", "Original"} <= set(agree.columns):
    a = agree.dropna()
    both = ((a.Toil < .05) & (a.Original < .05)).sum(); one = ((a.Toil < .05) ^ (a.Original < .05)).sum()
    print(f"\nSignificant in both pipelines: {both}; in only one pipeline: {one}; of {len(a)} comparable models")
print("Saved tables/Table3_survival_sensitivity.csv and tables/Table3_long.csv")

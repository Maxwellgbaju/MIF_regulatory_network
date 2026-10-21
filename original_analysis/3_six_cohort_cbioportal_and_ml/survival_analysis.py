import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test

surv_raw = pd.read_csv("survival_6cancers_raw.csv")
surv = surv_raw.pivot_table(index=["studyId","patientId"], columns="attribute", values="value", aggfunc="first").reset_index()

surv["event"] = surv["OS_STATUS"].astype(str).str.startswith("1").astype(int)
surv["time"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")
surv = surv.dropna(subset=["time"])

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
expr_wide = expr.pivot_table(index=["studyId","sampleId","patientId"], columns="geneSymbol", values="value").reset_index()

merged = expr_wide.merge(surv, on=["studyId","patientId"], how="inner")
print(f"Merged survival+expression shape: {merged.shape}\n")
print(merged[["studyId","time","event"]].groupby("studyId").agg(n=("time","count"), events=("event","sum")))
print()

genes = ["MIF", "SP1", "HIF1A", "RELA"]
results = []

for study in merged["studyId"].unique():
    sub = merged[merged["studyId"] == study]
    for gene in genes:
        d = sub[["time","event",gene]].dropna()
        if len(d) < 20 or d["event"].sum() < 5:
            print(f"Skipping {study}/{gene}: insufficient data (n={len(d)}, events={d['event'].sum()})")
            continue

        median_val = d[gene].median()
        high = d[d[gene] >= median_val]
        low = d[d[gene] < median_val]
        lr = logrank_test(high["time"], low["time"], high["event"], low["event"])

        cph = CoxPHFitter()
        cph.fit(d.rename(columns={gene: "expr"})[["time","event","expr"]], duration_col="time", event_col="event")
        hr = cph.hazard_ratios_["expr"]
        ci = cph.confidence_intervals_.loc["expr"]
        hr_ci_lo, hr_ci_hi = np.exp(ci.iloc[0]), np.exp(ci.iloc[1])
        cox_p = cph.summary.loc["expr", "p"]

        results.append({
            "study": study, "gene": gene, "n": len(d), "events": int(d["event"].sum()),
            "logrank_p": lr.p_value, "cox_HR": hr, "cox_HR_ci_lo": hr_ci_lo, "cox_HR_ci_hi": hr_ci_hi, "cox_p": cox_p
        })

res_df = pd.DataFrame(results)
res_df.to_csv("MIF_axis_survival_results.csv", index=False)
print("\n=== Survival results ===")
print(res_df.to_string(index=False))

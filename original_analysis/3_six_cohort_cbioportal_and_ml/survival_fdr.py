import pandas as pd
from statsmodels.stats.multitest import multipletests

res = pd.read_csv("MIF_axis_survival_results.csv")

# Primary test: Cox (continuous expression, no arbitrary cutpoint) -- FDR across all 24 tests
res["cox_fdr"] = multipletests(res["cox_p"], method="fdr_bh")[1]

# Secondary/consistency check: logrank (median-split KM) -- FDR across its own 24 tests
res["logrank_fdr"] = multipletests(res["logrank_p"], method="fdr_bh")[1]

# Flag agreement between the two tests at nominal p<0.05 (before FDR) -- a robustness signal
res["both_nominal_sig"] = (res["cox_p"] < 0.05) & (res["logrank_p"] < 0.05)

res = res.sort_values("cox_p")
res.to_csv("MIF_axis_survival_results_FDR.csv", index=False)

print(res[["study","gene","n","events","cox_HR","cox_p","cox_fdr","logrank_p","logrank_fdr","both_nominal_sig"]].to_string(index=False))

print("\n=== Surviving FDR < 0.05 (Cox) ===")
print(res[res["cox_fdr"] < 0.05][["study","gene","cox_HR","cox_p","cox_fdr"]].to_string(index=False))

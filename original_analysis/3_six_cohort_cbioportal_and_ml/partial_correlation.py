import pandas as pd
import numpy as np
from scipy.stats import spearmanr

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
proxy = pd.read_csv("purity_proxy_6cancers.csv")

expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
proxy_wide = proxy.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()

merged = expr_wide.merge(proxy_wide, on=["studyId","sampleId"], how="inner")
print(f"Merged shape: {merged.shape}")

def partial_spearman(x, y, z1, z2):
    # Rank-transform then compute linear residuals (Spearman partial correlation via rank + OLS residuals)
    from scipy.stats import rankdata
    rx, ry, rz1, rz2 = rankdata(x), rankdata(y), rankdata(z1), rankdata(z2)
    Z = np.column_stack([rz1, rz2, np.ones(len(rz1))])
    beta_x, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    beta_y, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    resid_x = rx - Z @ beta_x
    resid_y = ry - Z @ beta_y
    r, p = spearmanr(resid_x, resid_y)
    return r, p

results = []
for study in merged["studyId"].unique():
    sub = merged[merged["studyId"] == study].dropna(subset=["MIF","SP1","PTPRC","COL1A1"])
    if len(sub) < 10:
        continue
    raw_r, raw_p = spearmanr(sub["MIF"], sub["SP1"])
    partial_r, partial_p = partial_spearman(sub["MIF"].values, sub["SP1"].values, sub["PTPRC"].values, sub["COL1A1"].values)
    results.append({"study": study, "n": len(sub), "raw_rho": raw_r, "raw_p": raw_p,
                     "partial_rho": partial_r, "partial_p": partial_p})

res_df = pd.DataFrame(results)
res_df.to_csv("SP1_MIF_partial_correlation.csv", index=False)
print(res_df.to_string(index=False))

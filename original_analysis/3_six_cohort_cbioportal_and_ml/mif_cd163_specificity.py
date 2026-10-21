import pandas as pd
import numpy as np
from scipy.stats import spearmanr, rankdata

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
markers = pd.read_csv("immune_markers_6cancers.csv")
purity = pd.read_csv("purity_proxy_6cancers.csv")  # already has PTPRC from Phase 3.2

expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
markers_wide = markers.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
purity_wide = purity.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()

merged = expr_wide.merge(markers_wide, on=["studyId","sampleId"]).merge(purity_wide, on=["studyId","sampleId"])
print(f"Merged shape: {merged.shape}")

def partial_spearman(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    Z = np.column_stack([rz, np.ones(len(rz))])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    r, _ = spearmanr(rx - Z @ bx, ry - Z @ by)
    return r

results = []
for study in merged["studyId"].unique():
    sub = merged[["MIF", "CD163", "PTPRC"]].loc[merged["studyId"] == study].dropna()
    if len(sub) < 10:
        continue
    raw_r, raw_p = spearmanr(sub["MIF"], sub["CD163"])
    partial_r = partial_spearman(sub["MIF"].values, sub["CD163"].values, sub["PTPRC"].values)
    results.append({"study": study, "n": len(sub), "raw_rho": raw_r, "raw_p": raw_p, "partial_rho_ctrl_PTPRC": partial_r})

res_df = pd.DataFrame(results)
res_df.to_csv("MIF_CD163_specificity_check.csv", index=False)
print(res_df.to_string(index=False))

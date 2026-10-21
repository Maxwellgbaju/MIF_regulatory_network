import pandas as pd
import numpy as np
from scipy.stats import spearmanr, rankdata
from statsmodels.stats.multitest import multipletests

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
markers = pd.read_csv("immune_markers_6cancers.csv")
purity = pd.read_csv("purity_proxy_6cancers.csv")

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
    r, p = spearmanr(rx - Z @ bx, ry - Z @ by)
    return r, p

axis_genes = ["MIF", "SP1", "HIF1A", "RELA"]
immune_markers = ["PDCD1", "FOXP3", "CD163"]

results = []
for study in merged["studyId"].unique():
    for gene in axis_genes:
        for marker in immune_markers:
            sub = merged[[gene, marker, "PTPRC"]].loc[merged["studyId"] == study].dropna()
            if len(sub) < 10:
                continue
            raw_r, raw_p = spearmanr(sub[gene], sub[marker])
            partial_r, partial_p = partial_spearman(sub[gene].values, sub[marker].values, sub["PTPRC"].values)
            results.append({
                "study": study, "axis_gene": gene, "immune_marker": marker, "n": len(sub),
                "raw_rho": raw_r, "raw_p": raw_p,
                "partial_rho_ctrl_PTPRC": partial_r, "partial_p": partial_p
            })

res_df = pd.DataFrame(results)
res_df["raw_fdr"] = multipletests(res_df["raw_p"], method="fdr_bh")[1]
res_df["partial_fdr"] = multipletests(res_df["partial_p"], method="fdr_bh")[1]
res_df["attenuation_pct"] = 100 * (1 - res_df["partial_rho_ctrl_PTPRC"].abs() / res_df["raw_rho"].abs())
res_df = res_df.sort_values("partial_p")
res_df.to_csv("MIF_axis_immune_correlations_partial.csv", index=False)

print(f"\nTotal tests: {len(res_df)}")
print(f"Raw FDR<0.05: {(res_df['raw_fdr']<0.05).sum()}")
print(f"Partial FDR<0.05 (survives PTPRC adjustment): {(res_df['partial_fdr']<0.05).sum()}")

print("\n=== HIF1A - CD163 across all 6 cancers ===")
print(res_df[(res_df["axis_gene"]=="HIF1A") & (res_df["immune_marker"]=="CD163")]
      [["study","n","raw_rho","raw_fdr","partial_rho_ctrl_PTPRC","partial_fdr"]].to_string(index=False))

print("\n=== Full results, sorted by partial p-value ===")
print(res_df[["study","axis_gene","immune_marker","n","raw_rho","partial_rho_ctrl_PTPRC","partial_fdr"]].to_string(index=False))

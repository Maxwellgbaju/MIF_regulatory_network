import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
markers = pd.read_csv("immune_markers_6cancers.csv")

expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
markers_wide = markers.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
merged = expr_wide.merge(markers_wide, on=["studyId","sampleId"], how="inner")
print(f"Merged shape: {merged.shape}")

axis_genes = ["MIF", "SP1", "HIF1A", "RELA"]
immune_markers = ["PDCD1", "FOXP3", "CD163"]

results = []
for study in merged["studyId"].unique():
    sub = merged[merged["studyId"] == study]
    for gene in axis_genes:
        for marker in immune_markers:
            valid = sub[[gene, marker]].dropna()
            if len(valid) < 10:
                continue
            rho, pval = spearmanr(valid[gene], valid[marker])
            results.append({"study": study, "axis_gene": gene, "immune_marker": marker,
                             "n": len(valid), "rho": rho, "pval": pval})

res_df = pd.DataFrame(results)
res_df["fdr"] = multipletests(res_df["pval"], method="fdr_bh")[1]
res_df = res_df.sort_values("pval")
res_df.to_csv("MIF_axis_immune_correlations.csv", index=False)

print(f"\nTotal tests: {len(res_df)}")
print(res_df.to_string(index=False))

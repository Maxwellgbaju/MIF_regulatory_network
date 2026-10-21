import pandas as pd
from scipy.stats import spearmanr

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
markers = pd.read_csv("immune_markers_6cancers.csv")
purity = pd.read_csv("purity_proxy_6cancers.csv")

expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
markers_wide = markers.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
purity_wide = purity.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
merged = expr_wide.merge(markers_wide, on=["studyId","sampleId"]).merge(purity_wide, on=["studyId","sampleId"])

axis_genes = ["MIF", "SP1", "HIF1A", "RELA"]
immune_markers = ["PDCD1", "FOXP3", "CD163"]

print("=== Axis gene vs. overall immune content (PTPRC) ===")
for study in merged["studyId"].unique():
    for gene in axis_genes:
        sub = merged[[gene, "PTPRC"]].loc[merged["studyId"] == study].dropna()
        r, p = spearmanr(sub[gene], sub["PTPRC"])
        print(f"{study:35s} {gene:6s} vs PTPRC: rho={r:+.3f}  p={p:.2e}")

print("\n=== Immune marker vs. overall immune content (PTPRC) ===")
for study in merged["studyId"].unique():
    for marker in immune_markers:
        sub = merged[[marker, "PTPRC"]].loc[merged["studyId"] == study].dropna()
        r, p = spearmanr(sub[marker], sub["PTPRC"])
        print(f"{study:35s} {marker:6s} vs PTPRC: rho={r:+.3f}  p={p:.2e}")

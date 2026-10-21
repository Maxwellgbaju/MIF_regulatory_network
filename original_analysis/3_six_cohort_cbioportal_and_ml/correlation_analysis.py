import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

df = pd.read_csv("MIF_axis_expression_6cancers.csv")

# Pivot to wide format: one row per sample, one column per gene
wide = df.pivot_table(index=["studyId", "sampleId"], columns="geneSymbol", values="value").reset_index()
print(f"Wide-format shape: {wide.shape}")
print(wide.head())

studies = wide["studyId"].unique()
axis_genes = ["SP1", "HIF1A", "RELA"]

results = []
for study in studies:
    sub = wide[wide["studyId"] == study]
    for gene in axis_genes:
        valid = sub[["MIF", gene]].dropna()
        if len(valid) < 10:
            print(f"Skipping {study} / {gene}: only {len(valid)} paired samples")
            continue
        rho, pval = spearmanr(valid["MIF"], valid[gene])
        results.append({"study": study, "gene": gene, "n": len(valid), "rho": rho, "pval": pval})

res_df = pd.DataFrame(results)

# Benjamini-Hochberg FDR correction across ALL tests run (18 tests: 6 cancers x 3 genes)
res_df["fdr"] = multipletests(res_df["pval"], method="fdr_bh")[1]
res_df = res_df.sort_values("pval")

res_df.to_csv("MIF_axis_correlations_6cancers.csv", index=False)
print(f"\nSaved {len(res_df)} correlation tests to MIF_axis_correlations_6cancers.csv\n")
print(res_df.to_string(index=False))

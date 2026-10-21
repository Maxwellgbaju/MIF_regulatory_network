import json
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import pandas as pd

with open('xena_phenotype_full.json') as f:
    metadata = json.load(f)
with open('all_gene_expression.json') as f:
    expr_data = json.load(f)

genes = list(expr_data.keys())
print(f"Genes: {genes}")

# Identify all TCGA primary tumor cancer types
tcga_tumor_samples = {}  # disease -> [sample_ids]
for s, m in metadata.items():
    if m['_study'] == 'TCGA' and m['_sample_type'] == 'Primary Tumor':
        disease = m['primary disease or tissue']
        tcga_tumor_samples.setdefault(disease, []).append(s)

print(f"\nFound {len(tcga_tumor_samples)} TCGA cancer types with tumor samples")

# For each disease, determine its primary_site (from its own tumor samples), then find matched normal
results = []
for disease, tumor_ids in tcga_tumor_samples.items():
    if len(tumor_ids) < 10:
        continue  # skip cancer types with too few tumor samples

    sites = set(metadata[s]['_primary_site'] for s in tumor_ids if metadata[s]['_primary_site'])
    if not sites:
        continue
    primary_site = list(sites)[0]  # use the dominant/first site

    # Matched normal: TCGA Solid Tissue Normal (same disease) + GTEx Normal Tissue (same primary_site)
    tcga_normal_ids = [s for s, m in metadata.items()
                        if m['_study'] == 'TCGA' and m['_sample_type'] == 'Solid Tissue Normal'
                        and m['primary disease or tissue'] == disease]
    gtex_normal_ids = [s for s, m in metadata.items()
                        if m['_study'] == 'GTEX' and m['_sample_type'] == 'Normal Tissue'
                        and m['_primary_site'] == primary_site]
    normal_ids = tcga_normal_ids + gtex_normal_ids

    if len(normal_ids) < 10:
        continue  # skip if insufficient normal comparison

    for gene in genes:
        tumor_vals = np.array([float(expr_data[gene][s]) for s in tumor_ids if s in expr_data[gene]])
        normal_vals = np.array([float(expr_data[gene][s]) for s in normal_ids if s in expr_data[gene]])

        stat, pval = mannwhitneyu(tumor_vals, normal_vals, alternative='two-sided')
        log2fc = np.median(tumor_vals) - np.median(normal_vals)  # already log2 scale

        results.append({
            "gene": gene, "cancer_type": disease, "primary_site": primary_site,
            "n_tumor": len(tumor_vals), "n_normal_tcga": len(tcga_normal_ids), "n_normal_gtex": len(gtex_normal_ids),
            "median_tumor": np.median(tumor_vals), "median_normal": np.median(normal_vals),
            "log2FC": log2fc, "pval": pval
        })

res_df = pd.DataFrame(results)
res_df["FDR"] = multipletests(res_df["pval"], method="fdr_bh")[1]
res_df = res_df.sort_values("FDR")
res_df.to_csv("tumor_normal_full_results.csv", index=False)

print(f"\nTotal tests run: {len(res_df)}")
print(f"Cancer types with sufficient data: {res_df['cancer_type'].nunique()}")
print(f"\n=== Top 20 most significant gene-cancer results ===")
print(res_df.head(20).to_string(index=False))

print(f"\n=== MIF specifically, across all cancer types, ranked ===")
mif_results = res_df[res_df["gene"] == "MIF"].sort_values("FDR")
print(mif_results.to_string(index=False))

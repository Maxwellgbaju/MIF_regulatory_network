import json
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import pandas as pd

with open('lncrna_expression.json') as f:
    expr_data = json.load(f)
with open('xena_phenotype_full.json') as f:
    metadata = json.load(f)

gene = 'AL356488.2'
floor = -9.966

tcga_tumor_samples = {}
for s, m in metadata.items():
    if m['_study'] == 'TCGA' and m['_sample_type'] == 'Primary Tumor':
        disease = m['primary disease or tissue']
        tcga_tumor_samples.setdefault(disease, []).append(s)

results = []
for disease, tumor_ids in tcga_tumor_samples.items():
    if len(tumor_ids) < 10:
        continue
    sites = set(metadata[s]['_primary_site'] for s in tumor_ids if metadata[s]['_primary_site'])
    if not sites:
        continue
    primary_site = list(sites)[0]

    tcga_normal_ids = [s for s, m in metadata.items()
                        if m['_study'] == 'TCGA' and m['_sample_type'] == 'Solid Tissue Normal'
                        and m['primary disease or tissue'] == disease]
    gtex_normal_ids = [s for s, m in metadata.items()
                        if m['_study'] == 'GTEX' and m['_sample_type'] == 'Normal Tissue'
                        and m['_primary_site'] == primary_site]
    normal_ids = tcga_normal_ids + gtex_normal_ids
    if len(normal_ids) < 10:
        continue

    tumor_vals = [float(expr_data[gene][s]) for s in tumor_ids if s in expr_data[gene]]
    normal_vals = [float(expr_data[gene][s]) for s in normal_ids if s in expr_data[gene]]

    det_tumor = sum(1 for v in tumor_vals if v > floor)
    det_normal = sum(1 for v in normal_vals if v > floor)
    undet_tumor = len(tumor_vals) - det_tumor
    undet_normal = len(normal_vals) - det_normal

    table = [[det_tumor, undet_tumor], [det_normal, undet_normal]]
    odds_ratio, pval = fisher_exact(table)

    results.append({
        "cancer_type": disease, "primary_site": primary_site,
        "n_tumor": len(tumor_vals), "pct_detected_tumor": 100*det_tumor/len(tumor_vals),
        "n_normal": len(normal_vals), "pct_detected_normal": 100*det_normal/len(normal_vals),
        "odds_ratio": odds_ratio, "pval": pval
    })

res_df = pd.DataFrame(results)
res_df["FDR"] = multipletests(res_df["pval"], method="fdr_bh")[1]
res_df = res_df.sort_values("FDR")
res_df.to_csv("AL356488.2_detection_rate_results.csv", index=False)
print(res_df.to_string(index=False))

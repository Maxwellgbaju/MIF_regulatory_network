import xenaPython as xena
import json
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import pandas as pd

host = 'https://toil.xenahubs.net'
expr_dataset = 'TcgaTargetGtex_rsem_gene_tpm'

lncrna_ids = {
    "NORAD": "ENSG00000260032",
    "AL356488.2": "ENSG00000270066"
}

# Find exact versioned IDs (same verification step used for MIF/TFs)
all_ids = open('xena_all_gene_ids.txt').read().splitlines()
versioned_ids = {}
for symbol, base_id in lncrna_ids.items():
    matches = [x for x in all_ids if x.startswith(base_id)]
    if not matches:
        print(f"WARNING: {symbol} ({base_id}) not found in dataset -- skipping")
        continue
    versioned_ids[symbol] = matches[0]
    print(f"{symbol}: {matches[0]}")

with open('xena_phenotype_full.json') as f:
    metadata = json.load(f)
samples = list(metadata.keys())

print(f"\nFetching expression for {len(versioned_ids)} lncRNAs across {len(samples)} samples...")
result = xena.dataset_probe_values(host, expr_dataset, samples, list(versioned_ids.values()))
raw_values = result[1]

expr_data = {}
for idx, symbol in enumerate(versioned_ids.keys()):
    expr_data[symbol] = {samples[i]: raw_values[idx][i] for i in range(len(samples))}

with open("lncrna_expression.json", "w") as f:
    json.dump(expr_data, f)

for symbol in versioned_ids:
    vals = [float(v) for v in expr_data[symbol].values() if v not in (None, "NaN")]
    print(f"{symbol}: {len(vals)} valid values, range [{min(vals):.2f}, {max(vals):.2f}]")

# --- Same tumor-vs-normal grouping and statistics as the MIF/TF screen ---
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

    for gene in versioned_ids:
        tumor_vals = np.array([float(expr_data[gene][s]) for s in tumor_ids if s in expr_data[gene]])
        normal_vals = np.array([float(expr_data[gene][s]) for s in normal_ids if s in expr_data[gene]])
        stat, pval = mannwhitneyu(tumor_vals, normal_vals, alternative='two-sided')
        log2fc = np.median(tumor_vals) - np.median(normal_vals)
        results.append({
            "gene": gene, "cancer_type": disease, "primary_site": primary_site,
            "n_tumor": len(tumor_vals), "n_normal_tcga": len(tcga_normal_ids), "n_normal_gtex": len(gtex_normal_ids),
            "median_tumor": np.median(tumor_vals), "median_normal": np.median(normal_vals),
            "log2FC": log2fc, "pval": pval
        })

res_df = pd.DataFrame(results)
res_df["FDR"] = multipletests(res_df["pval"], method="fdr_bh")[1]
res_df = res_df.sort_values("FDR")
res_df.to_csv("lncrna_tumor_normal_results.csv", index=False)

print(f"\nTotal lncRNA tests run: {len(res_df)}")
print("\n=== Top 15 results ===")
print(res_df.head(15).to_string(index=False))

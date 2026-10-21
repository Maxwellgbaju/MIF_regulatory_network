import xenaPython as xena
import json
import numpy as np
import pandas as pd

host = 'https://toil.xenahubs.net'
expr_dataset = 'TcgaTargetGtex_rsem_gene_tpm'
surv_dataset = 'TCGA_survival_data'

with open('xena_phenotype_full.json') as f:
    metadata = json.load(f)
with open('all_gene_expression.json') as f:
    expr_data = json.load(f)
with open('lncrna_expression.json') as f:
    lncrna_data = json.load(f)

# Fetch RELA expression for all samples
samples = list(metadata.keys())
print("Fetching RELA expression...")
rela_values = xena.dataset_probe_values(host, expr_dataset, samples, ['ENSG00000173039.18'])
raw = rela_values[1][0]
rela_expr = {}
for i, s in enumerate(samples):
    try:
        v = float(raw[i])
    except (ValueError, TypeError):
        continue
    if np.isnan(v):
        continue
    rela_expr[s] = v
print(f"RELA: {len(rela_expr)} valid values")

# Combine all 15 candidate features
all_features = {**expr_data, **lncrna_data, "RELA": rela_expr}
genes = ["MIF", "RELA", "EGR1", "SP1", "NR1H2", "RXRA", "IRF1", "KLF4", "IRF2",
         "NFIL3", "ESR1", "TP53", "TFAP2A", "NORAD", "AL356488.2"]
print(f"Full candidate panel ({len(genes)} genes): {genes}")

# KIRC tumor samples
kirc_tumor_ids = [s for s, m in metadata.items()
                  if m['_study'] == 'TCGA' and m['_sample_type'] == 'Primary Tumor'
                  and m['primary disease or tissue'] == 'Kidney Clear Cell Carcinoma']
print(f"\nKIRC tumor samples: {len(kirc_tumor_ids)}")

# Fetch survival for KIRC tumor samples
surv_values = xena.dataset_probe_values(host, surv_dataset, kirc_tumor_ids, ['OS', 'OS.time'])
os_status_raw = surv_values[1][0]
os_time_raw = surv_values[1][1]

records = []
for i, s in enumerate(kirc_tumor_ids):
    try:
        os_stat = float(os_status_raw[i])
        os_days = float(os_time_raw[i])
    except (ValueError, TypeError):
        continue
    if np.isnan(os_stat) or np.isnan(os_days):
        continue

    row = {"sample": s, "OS": os_stat, "OS_time_months": os_days / 30.44}
    valid = True
    for gene in genes:
        if s not in all_features[gene]:
            valid = False
            break
        v = all_features[gene][s]
        try:
            v = float(v)
        except (ValueError, TypeError):
            valid = False
            break
        if np.isnan(v):
            valid = False
            break
        row[gene] = v
    if valid:
        records.append(row)

df = pd.DataFrame(records)
print(f"\nComplete KIRC dataset (all 15 genes + survival, no missing values): n={len(df)}, events={int(df['OS'].sum())}")
df.to_csv("capstone_kirc_full_matrix.csv", index=False)
print("Saved capstone_kirc_full_matrix.csv")
print(df[genes].describe())

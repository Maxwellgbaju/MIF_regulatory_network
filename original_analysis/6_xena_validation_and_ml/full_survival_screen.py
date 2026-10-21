import xenaPython as xena
import json
import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests

host = 'https://toil.xenahubs.net'
surv_dataset = 'TCGA_survival_data'

with open('xena_phenotype_full.json') as f:
    metadata = json.load(f)
with open('all_gene_expression.json') as f:
    expr_data = json.load(f)
with open('lncrna_expression.json') as f:
    lncrna_data = json.load(f)

all_expr = {**expr_data, **lncrna_data}
genes = list(all_expr.keys())
print(f"Genes to test: {genes}")

tcga_tumor_samples = {}
for s, m in metadata.items():
    if m['_study'] == 'TCGA' and m['_sample_type'] == 'Primary Tumor':
        disease = m['primary disease or tissue']
        tcga_tumor_samples.setdefault(disease, []).append(s)

all_tumor_ids = [s for ids in tcga_tumor_samples.values() for s in ids]
print(f"\nFetching survival data for {len(all_tumor_ids)} tumor samples...")
surv_values = xena.dataset_probe_values(host, surv_dataset, all_tumor_ids, ['OS', 'OS.time'])
os_status_raw = surv_values[1][0]
os_time_raw = surv_values[1][1]

survival = {}
n_rejected_nan = 0
for i, s in enumerate(all_tumor_ids):
    try:
        os_stat = float(os_status_raw[i])
        os_days = float(os_time_raw[i])
    except (ValueError, TypeError):
        continue
    # explicit NaN check: float() silently converts the string "NaN" without raising
    if np.isnan(os_stat) or np.isnan(os_days):
        n_rejected_nan += 1
        continue
    survival[s] = {"OS": os_stat, "OS_time_months": os_days / 30.44}

print(f"Samples with valid survival data: {len(survival)} (rejected {n_rejected_nan} with NaN values)")

results = []
for disease, tumor_ids in tcga_tumor_samples.items():
    valid_ids = [s for s in tumor_ids if s in survival]
    n_events = sum(survival[s]["OS"] for s in valid_ids)
    if len(valid_ids) < 20 or n_events < 5:
        continue

    for gene in genes:
        gene_vals = {}
        for s in valid_ids:
            if s not in all_expr[gene]:
                continue
            try:
                v = float(all_expr[gene][s])
            except (ValueError, TypeError):
                continue
            if np.isnan(v):
                continue
            gene_vals[s] = v

        if len(gene_vals) < 20:
            continue

        df = pd.DataFrame({
            "sample": list(gene_vals.keys()),
            "expr": list(gene_vals.values()),
        })
        df["OS"] = df["sample"].map(lambda s: survival[s]["OS"])
        df["OS_time"] = df["sample"].map(lambda s: survival[s]["OS_time_months"])

        median_val = df["expr"].median()
        high = df[df["expr"] >= median_val]
        low = df[df["expr"] < median_val]
        if len(high) < 10 or len(low) < 10:
            continue

        lr = logrank_test(high["OS_time"], low["OS_time"], high["OS"], low["OS"])

        try:
            cph = CoxPHFitter()
            cph.fit(df[["OS_time", "OS", "expr"]], duration_col="OS_time", event_col="OS")
            hr = cph.hazard_ratios_["expr"]
            ci = cph.confidence_intervals_.loc["expr"]
            hr_lo, hr_hi = np.exp(ci.iloc[0]), np.exp(ci.iloc[1])
        except Exception as e:
            print(f"COX FIT FAILED for {gene} in {disease}: {type(e).__name__}: {e}")
            hr, hr_lo, hr_hi = np.nan, np.nan, np.nan

        results.append({
            "gene": gene, "cancer_type": disease, "n_total": len(df), "n_events": int(df["OS"].sum()),
            "n_high": len(high), "n_low": len(low),
            "logrank_p": lr.p_value, "cox_HR": hr, "cox_HR_lo": hr_lo, "cox_HR_hi": hr_hi
        })

res_df = pd.DataFrame(results)
res_df["FDR"] = multipletests(res_df["logrank_p"], method="fdr_bh")[1]
res_df = res_df.sort_values("FDR")
res_df.to_csv("full_survival_screen_results.csv", index=False)

n_sig = (res_df["FDR"] < 0.05).sum()
n_nan_hr = res_df["cox_HR"].isna().sum()
print(f"\nTotal tests run: {len(res_df)}")
print(f"FDR-significant results (FDR < 0.05): {n_sig}")
print(f"Results with failed Cox fit (NaN HR): {n_nan_hr}")
print(f"\n=== Top 20 results by FDR ===")
print(res_df.head(20).to_string(index=False))

import pandas as pd
import numpy as np
from scipy.stats import spearmanr, rankdata

np.random.seed(42)
N_BOOT = 2000  # reduced from 5000 given 72x the workload; still gives stable estimates

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
markers = pd.read_csv("immune_markers_6cancers.csv")
purity = pd.read_csv("purity_proxy_6cancers.csv")

expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
markers_wide = markers.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
purity_wide = purity.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
merged = expr_wide.merge(markers_wide, on=["studyId","sampleId"]).merge(purity_wide, on=["studyId","sampleId"])

def partial_spearman(x, y, z):
    rx, ry, rz = rankdata(x), rankdata(y), rankdata(z)
    Z = np.column_stack([rz, np.ones(len(rz))])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    r, _ = spearmanr(rx - Z @ bx, ry - Z @ by)
    return r

axis_genes = ["MIF", "SP1", "HIF1A", "RELA"]
immune_markers = ["PDCD1", "FOXP3", "CD163"]

results = []
total = len(merged["studyId"].unique()) * len(axis_genes) * len(immune_markers)
count = 0
for study in merged["studyId"].unique():
    for gene in axis_genes:
        for marker in immune_markers:
            count += 1
            sub = merged[[gene, marker, "PTPRC"]].loc[merged["studyId"] == study].dropna().reset_index(drop=True)
            n = len(sub)
            if n < 10:
                continue

            obs_raw, _ = spearmanr(sub[gene], sub[marker])
            obs_partial = partial_spearman(sub[gene].values, sub[marker].values, sub["PTPRC"].values)
            obs_diff = obs_raw - obs_partial

            boot_diffs = np.empty(N_BOOT)
            idx_all = np.arange(n)
            for b in range(N_BOOT):
                idx = np.random.choice(idx_all, size=n, replace=True)
                bs = sub.iloc[idx]
                r_raw, _ = spearmanr(bs[gene], bs[marker])
                r_partial = partial_spearman(bs[gene].values, bs[marker].values, bs["PTPRC"].values)
                boot_diffs[b] = r_raw - r_partial

            ci_lo, ci_hi = np.percentile(boot_diffs, [2.5, 97.5])
            if obs_diff >= 0:
                p_boot = 2 * np.mean(boot_diffs <= 0)
            else:
                p_boot = 2 * np.mean(boot_diffs >= 0)
            p_boot = min(p_boot, 1.0)

            results.append({
                "study": study, "axis_gene": gene, "immune_marker": marker, "n": n,
                "raw_rho": obs_raw, "partial_rho": obs_partial, "diff": obs_diff,
                "diff_ci_lo": ci_lo, "diff_ci_hi": ci_hi, "bootstrap_p": p_boot
            })
            print(f"[{count}/{total}] {study[:4]} {gene}-{marker}: diff={obs_diff:+.3f}, boot_p={p_boot:.4f}")

res_df = pd.DataFrame(results)
res_df.to_csv("immune_bootstrap_diff_all.csv", index=False)
print(f"\nSaved {len(res_df)} rows to immune_bootstrap_diff_all.csv")
print(f"Significant raw-vs-partial shifts (bootstrap p<0.05): {(res_df['bootstrap_p']<0.05).sum()} / {len(res_df)}")

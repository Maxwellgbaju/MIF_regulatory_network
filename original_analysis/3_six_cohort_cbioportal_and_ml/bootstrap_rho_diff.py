import pandas as pd
import numpy as np
from scipy.stats import spearmanr, rankdata

np.random.seed(42)
N_BOOT = 5000

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
proxy = pd.read_csv("purity_proxy_6cancers.csv")

expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
proxy_wide = proxy.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
merged = expr_wide.merge(proxy_wide, on=["studyId","sampleId"], how="inner")

def partial_spearman(x, y, z1, z2):
    rx, ry, rz1, rz2 = rankdata(x), rankdata(y), rankdata(z1), rankdata(z2)
    Z = np.column_stack([rz1, rz2, np.ones(len(rz1))])
    bx, *_ = np.linalg.lstsq(Z, rx, rcond=None)
    by, *_ = np.linalg.lstsq(Z, ry, rcond=None)
    resid_x, resid_y = rx - Z @ bx, ry - Z @ by
    r, _ = spearmanr(resid_x, resid_y)
    return r

results = []
for study in merged["studyId"].unique():
    sub = merged[merged["studyId"] == study].dropna(subset=["MIF","SP1","PTPRC","COL1A1"]).reset_index(drop=True)
    n = len(sub)
    if n < 10:
        continue

    obs_raw, _ = spearmanr(sub["MIF"], sub["SP1"])
    obs_partial = partial_spearman(sub["MIF"].values, sub["SP1"].values, sub["PTPRC"].values, sub["COL1A1"].values)
    obs_diff = obs_raw - obs_partial

    boot_diffs = np.empty(N_BOOT)
    idx_all = np.arange(n)
    for b in range(N_BOOT):
        idx = np.random.choice(idx_all, size=n, replace=True)
        bs = sub.iloc[idx]
        r_raw, _ = spearmanr(bs["MIF"], bs["SP1"])
        r_partial = partial_spearman(bs["MIF"].values, bs["SP1"].values, bs["PTPRC"].values, bs["COL1A1"].values)
        boot_diffs[b] = r_raw - r_partial

    ci_lo, ci_hi = np.percentile(boot_diffs, [2.5, 97.5])
    # Two-sided bootstrap p-value: proportion of bootstrap diffs on the opposite side of zero from the observed diff, doubled
    if obs_diff >= 0:
        p_boot = 2 * np.mean(boot_diffs <= 0)
    else:
        p_boot = 2 * np.mean(boot_diffs >= 0)
    p_boot = min(p_boot, 1.0)

    results.append({
        "study": study, "n": n, "raw_rho": obs_raw, "partial_rho": obs_partial,
        "diff": obs_diff, "diff_ci_lo": ci_lo, "diff_ci_hi": ci_hi, "bootstrap_p": p_boot
    })
    print(f"{study}: diff={obs_diff:.4f}, 95% CI=[{ci_lo:.4f}, {ci_hi:.4f}], bootstrap p={p_boot:.4f}")

res_df = pd.DataFrame(results)
res_df.to_csv("SP1_MIF_bootstrap_diff_test.csv", index=False)
print("\nSaved SP1_MIF_bootstrap_diff_test.csv")

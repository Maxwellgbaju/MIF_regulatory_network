import pandas as pd
import numpy as np
from scipy.stats import spearmanr

np.random.seed(42)
N_BOOT = 2000

cna = pd.read_csv("cna_6cancers.csv")
expr = pd.read_csv("MIF_axis_expression_6cancers.csv")

cna_wide = cna.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
merged = cna_wide.merge(expr_wide, on=["studyId","sampleId"], suffixes=("_cna","_expr"))

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}

def rho_pair(df):
    all_r, _ = spearmanr(df["MIF_expr"], df["RELA_expr"])
    dip = df[(df["MIF_cna"] == 0) & (df["RELA_cna"] == 0)]
    dip_r, _ = spearmanr(dip["MIF_expr"], dip["RELA_expr"]) if len(dip) >= 10 else (np.nan, np.nan)
    return all_r, dip_r

# Also check MIF-RELA CNA co-occurrence (same logic as the SP1 check)
print("=== MIF vs RELA CNA co-occurrence check ===")
for study, code in study_to_code.items():
    sub = merged[merged["studyId"] == study][["MIF_cna","RELA_cna"]].dropna()
    r, p = spearmanr(sub["MIF_cna"], sub["RELA_cna"])
    print(f"{code:6s}: rho={r:+.3f}  p={p:.2e}  (n={len(sub)})")

print("\n=== All-samples vs CNA-diploid-only MIF-RELA correlation, with bootstrap significance ===")
results = []
for study, code in study_to_code.items():
    sub = merged[merged["studyId"] == study][["MIF_expr","RELA_expr","MIF_cna","RELA_cna"]].dropna().reset_index(drop=True)
    n = len(sub)
    obs_all, obs_dip = rho_pair(sub)
    obs_diff = obs_all - obs_dip

    boot_diffs = []
    idx_all = np.arange(n)
    for b in range(N_BOOT):
        bs = sub.iloc[np.random.choice(idx_all, size=n, replace=True)]
        r_all, r_dip = rho_pair(bs)
        if not np.isnan(r_dip):
            boot_diffs.append(r_all - r_dip)
    boot_diffs = np.array(boot_diffs)
    ci_lo, ci_hi = np.percentile(boot_diffs, [2.5, 97.5])
    p_boot = min(2 * np.mean(boot_diffs <= 0), 1.0) if obs_diff >= 0 else min(2 * np.mean(boot_diffs >= 0), 1.0)

    n_dip = len(sub[(sub["MIF_cna"]==0)&(sub["RELA_cna"]==0)])
    results.append({"cancer": code, "n_all": n, "rho_all": obs_all, "n_dip": n_dip, "rho_dip": obs_dip,
                     "diff": obs_diff, "diff_ci_lo": ci_lo, "diff_ci_hi": ci_hi, "bootstrap_p": p_boot})
    print(f"{code:6s}: all={obs_all:+.3f} (n={n})  |  diploid-only={obs_dip:+.3f} (n={n_dip})  |  diff={obs_diff:+.3f}  boot_p={p_boot:.4f}")

res_df = pd.DataFrame(results)
res_df.to_csv("RELA_CNA_diploid_bootstrap.csv", index=False)
print("\nSaved RELA_CNA_diploid_bootstrap.csv")

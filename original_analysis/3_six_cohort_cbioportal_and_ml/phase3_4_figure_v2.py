import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

np.random.seed(42)
N_BOOT = 2000

cna_freq = pd.read_csv("CNA_frequency_summary.csv")
mut_freq = pd.read_csv("Mutation_frequency_summary.csv")
cna = pd.read_csv("cna_6cancers.csv")
expr = pd.read_csv("MIF_axis_expression_6cancers.csv")

cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["MIF", "SP1", "HIF1A", "RELA"]
cancer_colors = dict(zip(cancer_order, ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b"]))

def fisher_ci(rho, n):
    z = np.arctanh(rho)
    se = 1 / np.sqrt(n - 3)
    zc = 1.959963984540054
    return np.tanh(z - zc*se), np.tanh(z + zc*se)

fig = plt.figure(figsize=(17, 9))
gs = fig.add_gridspec(2, 4, height_ratios=[1, 1], hspace=0.5, wspace=0.35)

# ---- Panel A: diverging bar chart, one small panel per gene, gain up / loss down ----
for idx, gene in enumerate(gene_order):
    ax = fig.add_subplot(gs[0, idx])
    sub = cna_freq[cna_freq["Gene"] == gene].set_index("Cancer").reindex(cancer_order)
    x = np.arange(len(cancer_order))
    ax.bar(x, sub["pct_any_gain"], color="#E69F00", label="% gain" if idx == 0 else None)
    ax.bar(x, -sub["pct_any_loss"], color="#5E4FA2", label="% loss" if idx == 0 else None)
    mut_sub = mut_freq[mut_freq["Gene"] == gene].set_index("Cancer").reindex(cancer_order)
    for xi, m in zip(x, mut_sub["pct_mutated"]):
        ax.text(xi, sub["pct_any_gain"].max() * 1.15 if sub["pct_any_gain"].max() > 0 else 5,
                f"{m:.1f}%", ha="center", fontsize=7, color="#333333")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(cancer_order, fontsize=8, rotation=45)
    ax.set_title(gene, fontsize=11)
    if idx == 0:
        ax.set_ylabel("% samples\n(gain above 0, loss below 0)", fontsize=9)
        ax.legend(fontsize=8, loc="upper right")

fig.text(0.13, 0.97, "A. Copy-number gain/loss frequency by gene (mutation % labeled above bars)", fontsize=12)

# ---- Panel B: all-samples vs diploid-only SP1-MIF correlation, with bootstrap significance ----
cna_wide = cna.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
expr_wide = expr.pivot_table(index=["studyId","sampleId"], columns="geneSymbol", values="value").reset_index()
merged = cna_wide.merge(expr_wide, on=["studyId","sampleId"], suffixes=("_cna","_expr"))
study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}

def rho_pair(df):
    all_r, _ = spearmanr(df["MIF_expr"], df["SP1_expr"])
    dip = df[(df["MIF_cna"] == 0) & (df["SP1_cna"] == 0)]
    dip_r, _ = spearmanr(dip["MIF_expr"], dip["SP1_expr"]) if len(dip) >= 10 else (np.nan, np.nan)
    return all_r, dip_r

results = []
for study, code in study_to_code.items():
    sub = merged[merged["studyId"] == study][["MIF_expr","SP1_expr","MIF_cna","SP1_cna"]].dropna().reset_index(drop=True)
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

    n_dip = len(sub[(sub["MIF_cna"]==0)&(sub["SP1_cna"]==0)])
    results.append({"cancer": code, "n_all": n, "rho_all": obs_all, "n_dip": n_dip, "rho_dip": obs_dip,
                     "diff": obs_diff, "diff_ci_lo": ci_lo, "diff_ci_hi": ci_hi, "bootstrap_p": p_boot})
    print(f"{code}: all={obs_all:.3f} dip={obs_dip:.3f} diff={obs_diff:+.3f} boot_p={p_boot:.4f}")

res_df = pd.DataFrame(results).set_index("cancer").reindex(cancer_order)
res_df.to_csv("CNA_diploid_bootstrap.csv")

ax1 = fig.add_subplot(gs[1, :])
x = np.arange(len(cancer_order)); width = 0.35
all_vals = res_df["rho_all"].values
dip_vals = res_df["rho_dip"].values
all_ci = [fisher_ci(res_df.loc[c,"rho_all"], res_df.loc[c,"n_all"]) for c in cancer_order]
dip_ci = [fisher_ci(res_df.loc[c,"rho_dip"], res_df.loc[c,"n_dip"]) for c in cancer_order]
all_err = [[all_vals[i]-all_ci[i][0] for i in range(6)],[all_ci[i][1]-all_vals[i] for i in range(6)]]
dip_err = [[dip_vals[i]-dip_ci[i][0] for i in range(6)],[dip_ci[i][1]-dip_vals[i] for i in range(6)]]

ax1.bar(x-width/2, all_vals, width, yerr=all_err, capsize=4, label="All samples", color="#C00000", error_kw={"elinewidth":1.2})
ax1.bar(x+width/2, dip_vals, width, yerr=dip_err, capsize=4, label="CNA-diploid only (MIF & SP1)", color="#2E7D32", error_kw={"elinewidth":1.2})

def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "n.s."

for i, c in enumerate(cancer_order):
    p_val = res_df.loc[c, "bootstrap_p"]
    lower = min(all_ci[i][0], dip_ci[i][0])
    bracket_y = lower - 0.05
    x1, x2 = x[i]-width/2, x[i]+width/2
    ax1.plot([x1,x1,x2,x2],[bracket_y+0.015,bracket_y,bracket_y,bracket_y+0.015], color="black", linewidth=1)
    ax1.text((x1+x2)/2, bracket_y-0.02, sig_label(p_val), ha="center", va="top", fontsize=9,
              fontstyle="italic" if sig_label(p_val)=="n.s." else "normal")

ax1.axhline(0, color="black", linewidth=0.8)
ax1.set_xticks(x); ax1.set_xticklabels(cancer_order, fontsize=10)
ax1.set_ylabel("Spearman ρ (SP1 vs MIF), 95% CI", fontsize=10)
ax1.set_ylim(bottom=min(all_ci[i][0] for i in range(6)) - 0.12, top=0.05)
ax1.set_title("B. SP1-MIF correlation persists without CNA (bootstrap test for shift significance)", fontsize=12, loc="left")
ax1.legend(fontsize=9, loc="lower right")

plt.savefig("Phase3_4_figure_v2.png", dpi=300, bbox_inches="tight")
plt.savefig("Phase3_4_figure_v2.eps", format="eps", bbox_inches="tight")
print("\nSaved Phase3_4_figure_v2.png/.eps")

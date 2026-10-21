import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

cna_freq = pd.read_csv("CNA_frequency_summary.csv")
mut_freq = pd.read_csv("Mutation_frequency_summary.csv")
sp1_diploid = pd.read_csv("CNA_diploid_bootstrap.csv")
rela_diploid = pd.read_csv("RELA_CNA_diploid_bootstrap.csv")

cancer_order = ["KIRC", "LIHC", "STAD", "BRCA", "LUSC", "SKCM"]
gene_order = ["MIF", "SP1", "HIF1A", "RELA"]

def fisher_ci(rho, n):
    z = np.arctanh(rho)
    se = 1 / np.sqrt(n - 3)
    zc = 1.959963984540054
    return np.tanh(z - zc*se), np.tanh(z + zc*se)

def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01: return "**"
    if p < 0.05: return "*"
    return "n.s."

fig = plt.figure(figsize=(17, 13))
gs = fig.add_gridspec(3, 4, height_ratios=[1, 1, 1], hspace=0.55, wspace=0.35)

# ---- Panel A: diverging bar chart, one small panel per gene ----
for idx, gene in enumerate(gene_order):
    ax = fig.add_subplot(gs[0, idx])
    sub = cna_freq[cna_freq["Gene"] == gene].set_index("Cancer").reindex(cancer_order)
    x = np.arange(len(cancer_order))
    ax.bar(x, sub["pct_any_gain"], color="#E69F00", label="% gain" if idx == 0 else None)
    ax.bar(x, -sub["pct_any_loss"], color="#5E4FA2", label="% loss" if idx == 0 else None)
    mut_sub = mut_freq[mut_freq["Gene"] == gene].set_index("Cancer").reindex(cancer_order)
    for xi, m in zip(x, mut_sub["pct_mutated"]):
        ax.text(xi, sub["pct_any_gain"].max() * 1.25 if sub["pct_any_gain"].max() > 0 else 5,
                f"{m:.1f}%", ha="center", fontsize=7, color="#333333")
    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(cancer_order, fontsize=8, rotation=45)
    ax.set_title(gene, fontsize=12, fontweight="bold", pad=18)
    if idx == 0:
        ax.set_ylabel("% samples\n(gain above 0, loss below 0)", fontsize=9)
        ax.legend(fontsize=7, loc="upper right", bbox_to_anchor=(1.0, 0.85))

fig.text(0.13, 0.90, "A. Copy-number gain/loss frequency by gene (mutation % labeled above bars)", fontsize=12)

def plot_diploid_panel(ax, df, gene_label, y_range=None):
    x = np.arange(len(cancer_order)); width = 0.35
    df_ordered = df.set_index("cancer").reindex(cancer_order)
    all_vals = df_ordered["rho_all"].values
    dip_vals = df_ordered["rho_dip"].values
    ns = df_ordered["n_all"].values
    n_dips = df_ordered["n_dip"].values

    all_ci = [fisher_ci(all_vals[i], ns[i]) for i in range(len(cancer_order))]
    dip_ci = [fisher_ci(dip_vals[i], n_dips[i]) for i in range(len(cancer_order))]
    all_err = [[all_vals[i]-all_ci[i][0] for i in range(6)],[all_ci[i][1]-all_vals[i] for i in range(6)]]
    dip_err = [[dip_vals[i]-dip_ci[i][0] for i in range(6)],[dip_ci[i][1]-dip_vals[i] for i in range(6)]]

    ax.bar(x-width/2, all_vals, width, yerr=all_err, capsize=4, label="All samples",
           color="#C00000", error_kw={"elinewidth":1.2})
    ax.bar(x+width/2, dip_vals, width, yerr=dip_err, capsize=4,
           label=f"CNA-diploid only\n(MIF & {gene_label})", color="#2E7D32", error_kw={"elinewidth":1.2})

    for i, c in enumerate(cancer_order):
        p_val = df_ordered.loc[c, "bootstrap_p"]
        lower = min(all_ci[i][0], dip_ci[i][0])
        bracket_y = lower - 0.06 if lower < 0 else -0.06
        x1, x2 = x[i]-width/2, x[i]+width/2
        ax.plot([x1,x1,x2,x2],[bracket_y+0.015,bracket_y,bracket_y,bracket_y+0.015], color="black", linewidth=1)
        label = sig_label(p_val)
        ax.text((x1+x2)/2, bracket_y-0.025, label, ha="center", va="top", fontsize=8,
                  fontstyle="italic" if label=="n.s." else "normal")

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x); ax.set_xticklabels(cancer_order, fontsize=9)
    ax.set_ylabel(f"Spearman ρ ({gene_label} vs MIF), 95% CI", fontsize=9)
    ax.legend(fontsize=8, loc="lower right" if all_vals.mean() < 0 else "upper right")

    if y_range is not None:
        ax.set_ylim(*y_range)
    else:
        y_lo = min(min(a[0] for a in all_ci), min(d[0] for d in dip_ci)) - 0.14
        y_hi = max(max(a[1] for a in all_ci), max(d[1] for d in dip_ci)) + 0.06
        ax.set_ylim(y_lo, y_hi)

# Shared symmetric y-range so the zero-line sits at the same height in both panels,
# making the opposite-sign relationship visually comparable rather than mismatched
shared_range = (-0.75, 0.75)

ax_b = fig.add_subplot(gs[1, :])
plot_diploid_panel(ax_b, sp1_diploid, "SP1", y_range=shared_range)
ax_b.set_title("B. SP1-MIF correlation persists without CNA (bootstrap test for shift significance)", fontsize=12, loc="left")

ax_c = fig.add_subplot(gs[2, :])
plot_diploid_panel(ax_c, rela_diploid, "RELA", y_range=shared_range)
ax_c.set_title("C. RELA-MIF correlation persists (and strengthens) without CNA (bootstrap test for shift significance)", fontsize=12, loc="left")

plt.savefig("Phase3_4_figure_final.png", dpi=300, bbox_inches="tight")
plt.savefig("Phase3_4_figure_final.eps", format="eps", bbox_inches="tight")
print("Saved Phase3_4_figure_final.png/.eps")

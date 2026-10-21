import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv
from lifelines.utils import concordance_index
from scipy.stats import wilcoxon

np.random.seed(42)
N_REPEATS = 20

# ============================================================
# STEP 1: Load data
# ============================================================
df = pd.read_csv("phase4_kirc_dataset.csv")
genes = ["MIF", "SP1", "HIF1A", "RELA"]
X = df[genes].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
y = Surv.from_dataframe("event", "time", df)
print(f"Dataset: n={len(df)}, events={df['event'].sum()}\n")

# ============================================================
# STEP 2: Select LASSO penalty via cross-validation
# ============================================================
print("=== Selecting LASSO penalty (alpha) via 5-fold CV ===")
alphas = np.logspace(-3, 0, 30)
kf_select = KFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = []
for alpha in alphas:
    fold_scores = []
    for train_idx, test_idx in kf_select.split(X_scaled):
        try:
            model = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[alpha], fit_baseline_model=True)
            model.fit(X_scaled[train_idx], y[train_idx])
            risk_test = model.predict(X_scaled[test_idx])
            c_idx = concordance_index(df["time"].values[test_idx], -risk_test, df["event"].values[test_idx])
            fold_scores.append(c_idx)
        except Exception:
            continue
    if fold_scores:
        cv_scores.append((alpha, np.mean(fold_scores)))
best_alpha, best_cv_score = max(cv_scores, key=lambda x: x[1])
print(f"Best alpha: {best_alpha:.5f}, CV C-index at selection: {best_cv_score:.4f}\n")

# ============================================================
# STEP 3: Refit on full data to get final coefficients
# ============================================================
final_model = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
final_model.fit(X_scaled, y)
coefs = final_model.coef_.flatten()
print("=== Final LASSO-selected coefficients (standardized) ===")
for gene, coef in zip(genes, coefs):
    status = "SELECTED" if coef != 0 else "dropped"
    print(f"  {gene}: {coef:+.4f}  [{status}]")
selected_genes = [g for g, c in zip(genes, coefs) if c != 0]
print(f"\nGenes retained by LASSO: {selected_genes}\n")

# ============================================================
# STEP 4: Repeated cross-validation, multi-gene score vs. MIF alone
# ============================================================
print(f"=== Repeated cross-validation ({N_REPEATS} x 5-fold = {N_REPEATS*5} paired estimates) ===")
multi_scores, mif_scores = [], []
for repeat in range(N_REPEATS):
    kf = KFold(n_splits=5, shuffle=True, random_state=repeat)
    for train_idx, test_idx in kf.split(X_scaled):
        m_multi = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
        m_multi.fit(X_scaled[train_idx], y[train_idx])
        risk_multi = m_multi.predict(X_scaled[test_idx])
        multi_scores.append(concordance_index(df["time"].values[test_idx], -risk_multi, df["event"].values[test_idx]))

        mif_train = X_scaled[train_idx, genes.index("MIF")].reshape(-1, 1)
        mif_test = X_scaled[test_idx, genes.index("MIF")].reshape(-1, 1)
        m_mif = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[1e-6], fit_baseline_model=True)
        m_mif.fit(mif_train, y[train_idx])
        risk_mif = m_mif.predict(mif_test)
        mif_scores.append(concordance_index(df["time"].values[test_idx], -risk_mif, df["event"].values[test_idx]))

multi_scores = np.array(multi_scores)
mif_scores = np.array(mif_scores)
diffs = multi_scores - mif_scores

print(f"LASSO-selected score: mean C-index = {multi_scores.mean():.4f} (sd={multi_scores.std():.4f})")
print(f"MIF alone:            mean C-index = {mif_scores.mean():.4f} (sd={mif_scores.std():.4f})")
print(f"Mean paired difference: {diffs.mean():+.4f}")

stat, wilcoxon_p = wilcoxon(diffs)
pct_favor_multi = (diffs > 0).mean() * 100
print(f"Wilcoxon signed-rank p-value: {wilcoxon_p:.4f}")
print(f"Fold-repeats favoring LASSO-selected score: {pct_favor_multi:.1f}%")

pd.DataFrame({"multi_gene_cindex": multi_scores, "mif_alone_cindex": mif_scores}).to_csv(
    "phase4_repeated_cv_scores.csv", index=False)
pd.DataFrame({"gene": genes, "standardized_coefficient": coefs}).to_csv(
    "phase4_risk_score_coefficients.csv", index=False)

# ============================================================
# STEP 5: Build the figure -- three panels, none resembling a KM curve
# ============================================================
fig = plt.figure(figsize=(15, 5.5))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.1, 1.1], wspace=0.35)

# --- Panel A: LASSO feature selection bar chart ---
ax0 = fig.add_subplot(gs[0])
colors = ["#888888" if c == 0 else "#C00000" for c in coefs]
bars = ax0.bar(genes, coefs, color=colors, edgecolor="black", linewidth=0.8)
ax0.axhline(0, color="black", linewidth=0.8)
for i, (gene, coef) in enumerate(zip(genes, coefs)):
    label = f"{coef:.3f}" if coef != 0 else "dropped"
    va = "bottom" if coef >= 0 else "top"
    offset = 0.005 if coef >= 0 else -0.005
    ax0.text(i, coef + offset, label, ha="center", va=va, fontsize=9,
              fontweight="bold" if coef != 0 else "normal",
              fontstyle="italic" if coef == 0 else "normal")
ax0.set_ylabel("Standardized LASSO-Cox coefficient")
ax0.set_title(f"A. Feature selection\n(alpha={best_alpha:.3f}, {len(selected_genes)}/{len(genes)} genes retained)", fontsize=10)

# --- Panel B: paired violin distribution of C-index across 100 fold-repeats ---
ax1 = fig.add_subplot(gs[1])
parts = ax1.violinplot([mif_scores, multi_scores], showmeans=True, showextrema=True, widths=0.7)
for pc, color in zip(parts['bodies'], ["#4472C4", "#C00000"]):
    pc.set_facecolor(color)
    pc.set_alpha(0.55)
    pc.set_edgecolor("black")
for partname in ('cbars', 'cmins', 'cmaxes', 'cmeans'):
    parts[partname].set_edgecolor("black")
ax1.set_xticks([1, 2])
ax1.set_xticklabels(["MIF alone", f"LASSO score\n({'+'.join(selected_genes)})"])
ax1.axhline(0.5, color="gray", linestyle=":", linewidth=1.2, label="Random (C=0.5)")
ax1.set_ylabel(f"C-index ({N_REPEATS*5} repeated-CV estimates)")
ax1.set_title("B. Cross-validated discrimination", fontsize=10)
ax1.legend(fontsize=8, loc="lower right")

# --- Panel C: paired difference histogram with Wilcoxon test annotated ---
ax2 = fig.add_subplot(gs[2])
ax2.hist(diffs, bins=18, color="#70AD47", edgecolor="black", alpha=0.85)
ax2.axvline(0, color="black", linestyle="--", linewidth=1.3)
ax2.axvline(diffs.mean(), color="#C00000", linestyle="-", linewidth=1.8,
            label=f"Mean = {diffs.mean():+.3f}")
stat_text = f"Wilcoxon p = {wilcoxon_p:.3f}\n{pct_favor_multi:.0f}% of folds favor LASSO score"
ax2.text(0.97, 0.95, stat_text, transform=ax2.transAxes, ha="right", va="top",
          fontsize=9, bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"))
ax2.set_xlabel("Paired C-index difference\n(LASSO score - MIF alone)")
ax2.set_ylabel("Count")
ax2.set_title("C. Paired difference per fold-repeat", fontsize=10)
ax2.legend(fontsize=8, loc="upper left")

plt.suptitle("Phase 4: LASSO-Cox risk score construction and cross-validated performance (KIRC)", fontsize=12, y=1.04)
plt.tight_layout()
plt.savefig("phase4_complete_figure.png", dpi=300, bbox_inches="tight")
plt.savefig("phase4_complete_figure.eps", format="eps", bbox_inches="tight")
print("\nSaved phase4_complete_figure.png AND phase4_complete_figure.eps")

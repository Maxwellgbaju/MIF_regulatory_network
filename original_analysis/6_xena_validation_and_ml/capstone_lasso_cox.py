import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv
from lifelines.utils import concordance_index
from scipy.stats import wilcoxon

np.random.seed(42)

df = pd.read_csv("capstone_kirc_full_matrix.csv")
genes = ["MIF", "RELA", "EGR1", "SP1", "NR1H2", "RXRA", "IRF1", "KLF4", "IRF2",
         "NFIL3", "ESR1", "TP53", "TFAP2A", "NORAD", "AL356488.2"]

X = df[genes].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
y = Surv.from_dataframe("OS", "OS_time_months", df)

print(f"Dataset: n={len(df)}, events={df['OS'].sum()}, features={len(genes)}\n")

# --- Step 1: select LASSO penalty via 5-fold CV ---
print("=== Selecting LASSO penalty (alpha) via 5-fold CV ===")
alphas = np.logspace(-3, 0, 40)
kf_select = KFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = []
for alpha in alphas:
    fold_scores = []
    for train_idx, test_idx in kf_select.split(X_scaled):
        try:
            model = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[alpha], fit_baseline_model=True)
            model.fit(X_scaled[train_idx], y[train_idx])
            risk_test = model.predict(X_scaled[test_idx])
            c_idx = concordance_index(df["OS_time_months"].values[test_idx], -risk_test, df["OS"].values[test_idx])
            fold_scores.append(c_idx)
        except Exception:
            continue
    if fold_scores:
        cv_scores.append((alpha, np.mean(fold_scores)))
best_alpha, best_cv_score = max(cv_scores, key=lambda x: x[1])
print(f"Best alpha: {best_alpha:.5f}, CV C-index at selection: {best_cv_score:.4f}\n")

# --- Step 2: refit on full data ---
final_model = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
final_model.fit(X_scaled, y)
coefs = final_model.coef_.flatten()

print("=== Final LASSO-selected coefficients (standardized) ===")
for gene, coef in sorted(zip(genes, coefs), key=lambda x: -abs(x[1])):
    status = "SELECTED" if coef != 0 else "dropped"
    print(f"  {gene:15s}: {coef:+.4f}  [{status}]")
selected_genes = [g for g, c in zip(genes, coefs) if c != 0]
print(f"\nGenes retained by LASSO ({len(selected_genes)}/{len(genes)}): {selected_genes}")

pd.DataFrame({"gene": genes, "standardized_coefficient": coefs}).to_csv(
    "capstone_lasso_coefficients.csv", index=False)

# --- Step 3: repeated CV, three-way comparison (multi-gene LASSO vs MIF alone vs RELA alone) ---
print(f"\n=== Repeated cross-validation (20 x 5-fold = 100 paired estimates) ===")
mif_idx = genes.index("MIF")
rela_idx = genes.index("RELA")

multi_scores, mif_scores, rela_scores = [], [], []
for repeat in range(20):
    kf = KFold(n_splits=5, shuffle=True, random_state=repeat)
    for train_idx, test_idx in kf.split(X_scaled):
        m_multi = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
        m_multi.fit(X_scaled[train_idx], y[train_idx])
        risk_multi = m_multi.predict(X_scaled[test_idx])
        multi_scores.append(concordance_index(df["OS_time_months"].values[test_idx], -risk_multi, df["OS"].values[test_idx]))

        for idx, score_list in [(mif_idx, mif_scores), (rela_idx, rela_scores)]:
            train_1g = X_scaled[train_idx, idx].reshape(-1, 1)
            test_1g = X_scaled[test_idx, idx].reshape(-1, 1)
            m_1g = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[1e-6], fit_baseline_model=True)
            m_1g.fit(train_1g, y[train_idx])
            risk_1g = m_1g.predict(test_1g)
            score_list.append(concordance_index(df["OS_time_months"].values[test_idx], -risk_1g, df["OS"].values[test_idx]))

multi_scores, mif_scores, rela_scores = map(np.array, [multi_scores, mif_scores, rela_scores])

print(f"Multi-gene LASSO score: mean C-index = {multi_scores.mean():.4f} (sd={multi_scores.std():.4f})")
print(f"MIF alone:              mean C-index = {mif_scores.mean():.4f} (sd={mif_scores.std():.4f})")
print(f"RELA alone:             mean C-index = {rela_scores.mean():.4f} (sd={rela_scores.std():.4f})")

for name_a, scores_a, name_b, scores_b in [
    ("Multi-gene", multi_scores, "MIF alone", mif_scores),
    ("Multi-gene", multi_scores, "RELA alone", rela_scores),
    ("RELA alone", rela_scores, "MIF alone", mif_scores),
]:
    diff = scores_a - scores_b
    stat, p = wilcoxon(diff)
    print(f"{name_a} vs {name_b}: mean diff={diff.mean():+.4f}, Wilcoxon p={p:.4f}, {name_a} wins {100*(diff>0).mean():.1f}% of folds")

pd.DataFrame({"multi_gene": multi_scores, "mif_alone": mif_scores, "rela_alone": rela_scores}).to_csv(
    "capstone_repeated_cv_scores.csv", index=False)
print("\nSaved capstone_lasso_coefficients.csv and capstone_repeated_cv_scores.csv")

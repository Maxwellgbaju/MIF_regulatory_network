import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import concordance_index

np.random.seed(42)

df = pd.read_csv("phase4_kirc_dataset.csv")
genes = ["MIF", "SP1", "HIF1A", "RELA"]

X = df[genes].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

y = Surv.from_dataframe("event", "time", df)

# --- Step 1: cross-validated selection of LASSO penalty strength (alpha) via C-index ---
print("=== Selecting LASSO penalty via 5-fold cross-validation ===")
alphas = np.logspace(-3, 0, 30)
kf = KFold(n_splits=5, shuffle=True, random_state=42)
cv_scores = []

for alpha in alphas:
    fold_scores = []
    for train_idx, test_idx in kf.split(X_scaled):
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
print(f"Best alpha: {best_alpha:.5f}, CV C-index: {best_cv_score:.4f}")

# --- Step 2: refit on full data with best alpha to get final coefficients ---
final_model = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
final_model.fit(X_scaled, y)
coefs = final_model.coef_.flatten()

print("\n=== Final risk-score coefficients (standardized) ===")
for gene, coef in zip(genes, coefs):
    print(f"  {gene}: {coef:+.4f}" + ("  (dropped by LASSO)" if coef == 0 else ""))

# --- Step 3: compute risk score for every patient, median split, KM curve ---
df["risk_score"] = X_scaled @ coefs
median_risk = df["risk_score"].median()
df["risk_group"] = np.where(df["risk_score"] >= median_risk, "High risk", "Low risk")

high = df[df["risk_group"] == "High risk"]
low = df[df["risk_group"] == "Low risk"]
lr = logrank_test(high["time"], low["time"], high["event"], low["event"])
print(f"\nRisk-score median split log-rank test: p = {lr.p_value:.2e}")

fig, ax = plt.subplots(figsize=(7, 6))
kmf_high = KaplanMeierFitter()
kmf_high.fit(high["time"], high["event"], label=f"High risk (n={len(high)})")
kmf_high.plot_survival_function(ax=ax, color="#C00000")
kmf_low = KaplanMeierFitter()
kmf_low.fit(low["time"], low["event"], label=f"Low risk (n={len(low)})")
kmf_low.plot_survival_function(ax=ax, color="#4472C4")
ax.set_xlabel("Time (months)")
ax.set_ylabel("Overall survival probability")
ax.set_title("KIRC: Overall survival by multi-gene risk score")
ax.text(0.98, 0.05, f"Log-rank p = {lr.p_value:.2e}", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=10, bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"))
ax.legend(loc="upper right")
plt.tight_layout()
plt.savefig("phase4_KM_riskscore.png", dpi=300, bbox_inches="tight")
print("Saved phase4_KM_riskscore.png")

# --- Step 4: cross-validated C-index comparison: multi-gene risk score vs MIF alone ---
print("\n=== Cross-validated C-index: multi-gene risk score vs. MIF alone ===")
cindex_multi, cindex_mif = [], []
for train_idx, test_idx in kf.split(X_scaled):
    # Multi-gene model (refit LASSO on this fold's training data at the selected alpha)
    m_multi = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
    m_multi.fit(X_scaled[train_idx], y[train_idx])
    risk_multi = m_multi.predict(X_scaled[test_idx])
    cindex_multi.append(concordance_index(df["time"].values[test_idx], -risk_multi, df["event"].values[test_idx]))

    # MIF-only model
    mif_train = X_scaled[train_idx, genes.index("MIF")].reshape(-1, 1)
    mif_test = X_scaled[test_idx, genes.index("MIF")].reshape(-1, 1)
    m_mif = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[1e-6], fit_baseline_model=True)
    m_mif.fit(mif_train, y[train_idx])
    risk_mif = m_mif.predict(mif_test)
    cindex_mif.append(concordance_index(df["time"].values[test_idx], -risk_mif, df["event"].values[test_idx]))

print(f"Multi-gene risk score: mean CV C-index = {np.mean(cindex_multi):.4f} (fold scores: {[f'{c:.3f}' for c in cindex_multi]})")
print(f"MIF alone:             mean CV C-index = {np.mean(cindex_mif):.4f} (fold scores: {[f'{c:.3f}' for c in cindex_mif]})")

results_summary = pd.DataFrame({
    "gene": genes, "standardized_coefficient": coefs
})
results_summary.to_csv("phase4_risk_score_coefficients.csv", index=False)
df.to_csv("phase4_kirc_with_risk_scores.csv", index=False)
print("\nSaved phase4_risk_score_coefficients.csv and phase4_kirc_with_risk_scores.csv")

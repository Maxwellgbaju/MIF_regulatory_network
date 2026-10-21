import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv
from sksurv.metrics import cumulative_dynamic_auc

np.random.seed(42)

df = pd.read_csv("phase4_kirc_dataset.csv")
genes = ["MIF", "SP1", "HIF1A", "RELA"]
X = df[genes].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
y = Surv.from_dataframe("event", "time", df)

best_alpha = 0.05736  # from Phase 4.1 selection

# Time points for evaluation: 12, 36, 60 months (1yr, 3yr, 5yr), within the observed follow-up range
max_time = df["time"].max()
eval_times = np.array([t for t in [12, 24, 36, 48, 60] if t < max_time])
print(f"Evaluating at time points (months): {eval_times}")

# --- 5-fold CV: compute time-dependent AUC(t) for both models on held-out folds ---
kf = KFold(n_splits=5, shuffle=True, random_state=42)
auc_multi_folds, auc_mif_folds = [], []

for fold_i, (train_idx, test_idx) in enumerate(kf.split(X_scaled)):
    y_train, y_test = y[train_idx], y[test_idx]

    # Multi-gene (LASSO-selected) model
    m_multi = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
    m_multi.fit(X_scaled[train_idx], y_train)
    risk_multi_test = m_multi.predict(X_scaled[test_idx])

    # MIF-only model
    mif_train = X_scaled[train_idx, genes.index("MIF")].reshape(-1, 1)
    mif_test = X_scaled[test_idx, genes.index("MIF")].reshape(-1, 1)
    m_mif = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[1e-6], fit_baseline_model=True)
    m_mif.fit(mif_train, y_train)
    risk_mif_test = m_mif.predict(mif_test)

    try:
        auc_multi, _ = cumulative_dynamic_auc(y_train, y_test, risk_multi_test, eval_times)
        auc_mif, _ = cumulative_dynamic_auc(y_train, y_test, risk_mif_test, eval_times)
        auc_multi_folds.append(auc_multi)
        auc_mif_folds.append(auc_mif)
    except Exception as e:
        print(f"Fold {fold_i}: skipped ({e})")

auc_multi_folds = np.array(auc_multi_folds)
auc_mif_folds = np.array(auc_mif_folds)

mean_auc_multi = auc_multi_folds.mean(axis=0)
mean_auc_mif = auc_mif_folds.mean(axis=0)
se_auc_multi = auc_multi_folds.std(axis=0) / np.sqrt(len(auc_multi_folds))
se_auc_mif = auc_mif_folds.std(axis=0) / np.sqrt(len(auc_mif_folds))

print("\n=== Time-dependent AUC (mean across 5 CV folds) ===")
for t, am, sm, ai, si in zip(eval_times, mean_auc_multi, se_auc_multi, mean_auc_mif, se_auc_mif):
    print(f"  t={t:.0f} months: LASSO score AUC={am:.3f} (SE={sm:.3f})  |  MIF alone AUC={ai:.3f} (SE={si:.3f})")

results_df = pd.DataFrame({
    "time_months": eval_times,
    "auc_lasso_score": mean_auc_multi, "se_lasso_score": se_auc_multi,
    "auc_mif_alone": mean_auc_mif, "se_mif_alone": se_auc_mif,
})
results_df.to_csv("phase4_2_time_dependent_auc.csv", index=False)

# --- Figure: AUC(t) curves for both models with SE bands ---
fig, ax = plt.subplots(figsize=(8, 6))
ax.plot(eval_times, mean_auc_multi, marker="o", color="#C00000", label="LASSO score (RELA)", linewidth=2)
ax.fill_between(eval_times, mean_auc_multi - se_auc_multi, mean_auc_multi + se_auc_multi, color="#C00000", alpha=0.2)
ax.plot(eval_times, mean_auc_mif, marker="s", color="#4472C4", label="MIF alone", linewidth=2)
ax.fill_between(eval_times, mean_auc_mif - se_auc_mif, mean_auc_mif + se_auc_mif, color="#4472C4", alpha=0.2)
ax.axhline(0.5, color="gray", linestyle=":", linewidth=1.2, label="Random (AUC=0.5)")
ax.set_xlabel("Time (months)")
ax.set_ylabel("Time-dependent AUC (5-fold CV mean ± SE)")
ax.set_title("KIRC: Time-dependent discrimination, LASSO score vs. MIF alone")
ax.legend(loc="best")
ax.set_ylim(0.4, 0.8)
plt.tight_layout()
plt.savefig("phase4_2_time_dependent_auc.png", dpi=300, bbox_inches="tight")
plt.savefig("phase4_2_time_dependent_auc.eps", format="eps", bbox_inches="tight")
print("\nSaved phase4_2_time_dependent_auc.png/.eps and phase4_2_time_dependent_auc.csv")

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
from scipy.stats import wilcoxon

np.random.seed(42)
N_REPEATS = 20

df = pd.read_csv("phase4_kirc_dataset.csv")
genes = ["MIF", "SP1", "HIF1A", "RELA"]
X = df[genes].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
y = Surv.from_dataframe("event", "time", df)

best_alpha = 0.05736
max_time = df["time"].max()
eval_times = np.array([t for t in [12, 24, 36, 48, 60] if t < max_time])
print(f"Evaluating at time points (months): {eval_times}")
print(f"Running {N_REPEATS} repeats x 5-fold CV = {N_REPEATS*5} estimates per time point per model\n")

all_auc_multi = []  # shape: (n_estimates, n_timepoints)
all_auc_mif = []

for repeat in range(N_REPEATS):
    kf = KFold(n_splits=5, shuffle=True, random_state=repeat)
    for train_idx, test_idx in kf.split(X_scaled):
        y_train, y_test = y[train_idx], y[test_idx]

        m_multi = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
        m_multi.fit(X_scaled[train_idx], y_train)
        risk_multi_test = m_multi.predict(X_scaled[test_idx])

        mif_train = X_scaled[train_idx, genes.index("MIF")].reshape(-1, 1)
        mif_test = X_scaled[test_idx, genes.index("MIF")].reshape(-1, 1)
        m_mif = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[1e-6], fit_baseline_model=True)
        m_mif.fit(mif_train, y_train)
        risk_mif_test = m_mif.predict(mif_test)

        try:
            auc_multi, _ = cumulative_dynamic_auc(y_train, y_test, risk_multi_test, eval_times)
            auc_mif, _ = cumulative_dynamic_auc(y_train, y_test, risk_mif_test, eval_times)
            all_auc_multi.append(auc_multi)
            all_auc_mif.append(auc_mif)
        except Exception:
            continue

all_auc_multi = np.array(all_auc_multi)
all_auc_mif = np.array(all_auc_mif)
print(f"Successful fold-repeats: {len(all_auc_multi)} / {N_REPEATS*5}\n")

results = []
print("=== Time-dependent AUC: repeated CV mean, SE, and paired significance ===")
for i, t in enumerate(eval_times):
    multi_t = all_auc_multi[:, i]
    mif_t = all_auc_mif[:, i]
    diff_t = multi_t - mif_t
    stat, p = wilcoxon(diff_t)
    pct_favor = (diff_t > 0).mean() * 100
    print(f"t={t:.0f}mo: LASSO={multi_t.mean():.3f}(SE={multi_t.std()/np.sqrt(len(multi_t)):.3f})  "
          f"MIF={mif_t.mean():.3f}(SE={mif_t.std()/np.sqrt(len(mif_t)):.3f})  "
          f"diff={diff_t.mean():+.3f}  wilcoxon_p={p:.4f}  %favor_LASSO={pct_favor:.0f}%")
    results.append({"time_months": t, "auc_lasso_mean": multi_t.mean(), "auc_lasso_se": multi_t.std()/np.sqrt(len(multi_t)),
                     "auc_mif_mean": mif_t.mean(), "auc_mif_se": mif_t.std()/np.sqrt(len(mif_t)),
                     "mean_diff": diff_t.mean(), "wilcoxon_p": p, "pct_favor_lasso": pct_favor})

results_df = pd.DataFrame(results)
results_df.to_csv("phase4_2_repeated_time_auc.csv", index=False)

fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(eval_times, results_df["auc_lasso_mean"], yerr=results_df["auc_lasso_se"],
            marker="o", color="#C00000", label="LASSO score (RELA)", linewidth=2, capsize=4)
ax.errorbar(eval_times, results_df["auc_mif_mean"], yerr=results_df["auc_mif_se"],
            marker="s", color="#4472C4", label="MIF alone", linewidth=2, capsize=4)
ax.axhline(0.5, color="gray", linestyle=":", linewidth=1.2, label="Random (AUC=0.5)")
for _, row in results_df.iterrows():
    sig = "*" if row["wilcoxon_p"] < 0.05 else ""
    if sig:
        ax.text(row["time_months"], max(row["auc_lasso_mean"], row["auc_mif_mean"]) + 0.02, sig,
                ha="center", fontsize=14, fontweight="bold")
ax.set_xlabel("Time (months)")
ax.set_ylabel(f"Time-dependent AUC ({N_REPEATS}x5-fold CV mean ± SE)")
ax.set_title("KIRC: Time-dependent discrimination, repeated CV")
ax.legend(loc="best")
ax.set_ylim(0.4, 0.8)
plt.tight_layout()
plt.savefig("phase4_2_repeated_time_auc.png", dpi=300, bbox_inches="tight")
plt.savefig("phase4_2_repeated_time_auc.eps", format="eps", bbox_inches="tight")
print("\nSaved phase4_2_repeated_time_auc.png/.eps and .csv")

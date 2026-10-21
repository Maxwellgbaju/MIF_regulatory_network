import pandas as pd
import numpy as np
from sksurv.ensemble import RandomSurvivalForest
from sksurv.util import Surv
from sksurv.metrics import as_concordance_index_ipcw_scorer
from sklearn.inspection import permutation_importance
from sklearn.model_selection import KFold
from lifelines.utils import concordance_index

np.random.seed(42)

df = pd.read_csv("capstone_cbioportal_kirc_full.csv")
genes = ["MIF", "RELA", "EGR1", "SP1", "NR1H2", "RXRA", "IRF1", "KLF4", "IRF2",
         "NFIL3", "ESR1", "TP53", "TFAP2A", "NORAD"]

X = df[genes].values
y = Surv.from_dataframe("event", "time", df)

print(f"Dataset: n={len(df)}, events={df['event'].sum()}\n")

rsf = RandomSurvivalForest(n_estimators=500, min_samples_split=10, min_samples_leaf=15,
                             max_features="sqrt", n_jobs=-1, random_state=42)
rsf.fit(X, y)

print("Computing permutation importance (this may take a minute)...")
result = permutation_importance(rsf, X, y, n_repeats=15, random_state=42, n_jobs=-1)

print("\n=== Random Survival Forest permutation importance ===")
for gene, imp, std in sorted(zip(genes, result.importances_mean, result.importances_std), key=lambda x: -x[1]):
    print(f"  {gene:15s}: {imp:+.4f} (sd={std:.4f})")

pd.DataFrame({"gene": genes, "importance_mean": result.importances_mean, "importance_std": result.importances_std}).sort_values(
    "importance_mean", ascending=False).to_csv("capstone_rsf_importance.csv", index=False)

print(f"\n=== Repeated cross-validation (10 x 5-fold = 50 estimates) ===")
rsf_scores = []
for repeat in range(10):
    kf = KFold(n_splits=5, shuffle=True, random_state=repeat)
    for train_idx, test_idx in kf.split(X):
        m = RandomSurvivalForest(n_estimators=300, min_samples_split=10, min_samples_leaf=15,
                                   max_features="sqrt", n_jobs=-1, random_state=repeat)
        m.fit(X[train_idx], y[train_idx])
        risk = m.predict(X[test_idx])
        rsf_scores.append(concordance_index(df["time"].values[test_idx], -risk, df["event"].values[test_idx]))

rsf_scores = np.array(rsf_scores)
print(f"RSF (all 14 genes): mean C-index = {rsf_scores.mean():.4f} (sd={rsf_scores.std():.4f})")
print(f"(For comparison: LASSO 9-gene score was 0.6143, MIF alone 0.5719, RELA alone 0.5901)")

pd.DataFrame({"rsf_cindex": rsf_scores}).to_csv("capstone_rsf_cv_scores.csv", index=False)
print("\nSaved capstone_rsf_importance.csv and capstone_rsf_cv_scores.csv")

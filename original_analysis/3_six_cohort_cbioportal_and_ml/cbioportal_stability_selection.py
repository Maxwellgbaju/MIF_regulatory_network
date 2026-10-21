import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv

np.random.seed(42)

df = pd.read_csv("capstone_cbioportal_kirc_full.csv")
genes = ["MIF", "RELA", "EGR1", "SP1", "NR1H2", "RXRA", "IRF1", "KLF4", "IRF2",
         "NFIL3", "ESR1", "TP53", "TFAP2A", "NORAD"]

X = df[genes].values
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)
y = Surv.from_dataframe("event", "time", df)
best_alpha = 0.01701  # from prior selection

selection_counts = {g: 0 for g in genes}
total_fits = 0
for repeat in range(20):
    kf = KFold(n_splits=5, shuffle=True, random_state=repeat)
    for train_idx, _ in kf.split(X_scaled):
        m = CoxnetSurvivalAnalysis(l1_ratio=1.0, alphas=[best_alpha], fit_baseline_model=True)
        m.fit(X_scaled[train_idx], y[train_idx])
        coefs = m.coef_.flatten()
        for g, c in zip(genes, coefs):
            if c != 0:
                selection_counts[g] += 1
        total_fits += 1

print(f"Total resampled fits: {total_fits}")
result = pd.DataFrame({"gene": genes,
                        "selection_pct": [100*selection_counts[g]/total_fits for g in genes]}).sort_values("selection_pct", ascending=False)
print(result.to_string(index=False))
result.to_csv("cbioportal_stability_selection.csv", index=False)

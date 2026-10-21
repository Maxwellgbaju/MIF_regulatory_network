import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from statsmodels.stats.multitest import multipletests

df = pd.read_csv("phase4_2_repeated_time_auc.csv")
df["fdr_p"] = multipletests(df["wilcoxon_p"], method="fdr_bh")[1]
df.to_csv("phase4_2_repeated_time_auc_FDR.csv", index=False)

print(df[["time_months","auc_lasso_mean","auc_mif_mean","mean_diff","wilcoxon_p","fdr_p"]].to_string(index=False))

fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(df["time_months"], df["auc_lasso_mean"], yerr=df["auc_lasso_mean"] - df["auc_lasso_mean"],
            marker="o", color="#C00000", label="LASSO score (RELA)", linewidth=2, capsize=4)
ax.errorbar(df["time_months"], df["auc_mif_mean"], marker="s", color="#4472C4", label="MIF alone", linewidth=2, capsize=4)

# Use standard errors from the original repeated CV (re-derive approx from data if not present)

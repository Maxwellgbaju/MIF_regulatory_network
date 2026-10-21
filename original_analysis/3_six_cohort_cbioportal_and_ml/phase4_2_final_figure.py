import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

df = pd.read_csv("phase4_2_repeated_time_auc_FDR.csv")

fig, ax = plt.subplots(figsize=(8, 6))
ax.errorbar(df["time_months"], df["auc_lasso_mean"], yerr=df["auc_lasso_se"],
            marker="o", color="#C00000", label="LASSO score (RELA)", linewidth=2, capsize=4)
ax.errorbar(df["time_months"], df["auc_mif_mean"], yerr=df["auc_mif_se"],
            marker="s", color="#4472C4", label="MIF alone", linewidth=2, capsize=4)
ax.axhline(0.5, color="gray", linestyle=":", linewidth=1.2, label="Random (AUC=0.5)")

for _, row in df.iterrows():
    if row["fdr_p"] < 0.05:
        y_pos = max(row["auc_lasso_mean"], row["auc_mif_mean"]) + 0.025
        ax.text(row["time_months"], y_pos, "*", ha="center", fontsize=16, fontweight="bold", color="black")
    ax.text(row["time_months"], 0.42, f"FDR p={row['fdr_p']:.3f}", ha="center", fontsize=7, rotation=90, color="#555555")

ax.set_xlabel("Time (months)")
ax.set_ylabel(f"Time-dependent AUC (20x5-fold CV mean ± SE)")
ax.set_title("KIRC: Time-dependent discrimination, RELA vs. MIF alone\n(* = FDR-significant, BH correction across 5 time points)")
ax.legend(loc="upper right")
ax.set_ylim(0.4, 0.75)
plt.tight_layout()
plt.savefig("phase4_2_final_figure.png", dpi=300, bbox_inches="tight")
plt.savefig("phase4_2_final_figure.eps", format="eps", bbox_inches="tight")
print("Saved phase4_2_final_figure.png/.eps")

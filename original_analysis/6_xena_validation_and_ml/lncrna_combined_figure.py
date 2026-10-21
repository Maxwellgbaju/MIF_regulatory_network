import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# --- NORAD data: filter from the full lncrna_tumor_normal_results.csv ---
norad_df = pd.read_csv("lncrna_tumor_normal_results.csv")
norad_df = norad_df[norad_df["gene"] == "NORAD"].sort_values("log2FC")

# --- AL356488.2 data: from the detection-rate analysis ---
al_df = pd.read_csv("AL356488.2_detection_rate_results.csv")
def corrected_log2_or(row):
    det_t = row["n_tumor"] * row["pct_detected_tumor"] / 100
    undet_t = row["n_tumor"] - det_t
    det_n = row["n_normal"] * row["pct_detected_normal"] / 100
    undet_n = row["n_normal"] - det_n
    if min(det_t, undet_t, det_n, undet_n) == 0:
        det_t, undet_t, det_n, undet_n = det_t+0.5, undet_t+0.5, det_n+0.5, undet_n+0.5
    return np.log2((det_t * undet_n) / (undet_t * det_n))
al_df["log2OR"] = al_df.apply(corrected_log2_or, axis=1)
al_df = al_df.sort_values("log2OR")

def sig_label(fdr):
    if fdr < 0.001: return "***"
    if fdr < 0.01: return "**"
    if fdr < 0.05: return "*"
    return ""

fig, axes = plt.subplots(1, 2, figsize=(17, 10))

# Panel A: NORAD
ax = axes[0]
colors = ["#C00000" if fdr < 0.05 else "#BBBBBB" for fdr in norad_df["FDR"]]
ax.barh(norad_df["cancer_type"], norad_df["log2FC"], color=colors, edgecolor="black", linewidth=0.6)
for i, (val, fdr) in enumerate(zip(norad_df["log2FC"], norad_df["FDR"])):
    label = sig_label(fdr)
    offset = 0.05 if val >= 0 else -0.05
    ha = "left" if val >= 0 else "right"
    ax.text(val + offset, i, label, va="center", ha=ha, fontsize=10, fontweight="bold")
ax.axvline(0, color="black", linewidth=1)
ax.set_xlabel("log2 fold-change (median expression, Tumor vs. Normal)", fontsize=10)
ax.set_title("A. NORAD", fontsize=13, fontweight="bold", loc="left")
ax.legend(handles=[Patch(facecolor="#C00000", edgecolor="black", label="FDR < 0.05"),
                    Patch(facecolor="#BBBBBB", edgecolor="black", label="Not significant")],
          loc="lower right", fontsize=9)

# Panel B: AL356488.2
ax2 = axes[1]
colors2 = ["#1565C0" if fdr < 0.05 else "#BBBBBB" for fdr in al_df["FDR"]]
ax2.barh(al_df["cancer_type"], al_df["log2OR"], color=colors2, edgecolor="black", linewidth=0.6)
for i, (val, fdr) in enumerate(zip(al_df["log2OR"], al_df["FDR"])):
    label = sig_label(fdr)
    offset = 0.15 if val >= 0 else -0.15
    ha = "left" if val >= 0 else "right"
    ax2.text(val + offset, i, label, va="center", ha=ha, fontsize=10, fontweight="bold")
ax2.axvline(0, color="black", linewidth=1)
ax2.set_xlabel("log2(Odds Ratio), continuity-corrected\n(% samples with detectable expression, Tumor vs. Normal)", fontsize=10)
ax2.set_title("B. AL356488.2", fontsize=13, fontweight="bold", loc="left")
ax2.legend(handles=[Patch(facecolor="#1565C0", edgecolor="black", label="FDR < 0.05"),
                     Patch(facecolor="#BBBBBB", edgecolor="black", label="Not significant")],
           loc="lower right", fontsize=9)

plt.suptitle("Pan-cancer tumor-vs-normal comparison for two candidate lncRNA sponges\n" +
             "(A) NORAD: continuous expression fold-change  |  (B) AL356488.2: detection-rate odds ratio (different statistic, see Methods)",
             fontsize=12, y=1.02)
plt.tight_layout()
plt.savefig("lncRNA_combined_figure.png", dpi=300, bbox_inches="tight")
plt.savefig("lncRNA_combined_figure.eps", format="eps", bbox_inches="tight")
print("Saved lncRNA_combined_figure.png/.eps")

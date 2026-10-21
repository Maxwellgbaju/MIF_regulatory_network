import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
from lifelines.plotting import add_at_risk_counts

# --- Rebuild the merged KIRC/RELA survival dataset ---
surv_raw = pd.read_csv("survival_6cancers_raw.csv")
surv = surv_raw.pivot_table(index=["studyId","patientId"], columns="attribute", values="value", aggfunc="first").reset_index()
surv["event"] = surv["OS_STATUS"].astype(str).str.startswith("1").astype(int)
surv["time"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")
surv = surv.dropna(subset=["time"])

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
expr_wide = expr.pivot_table(index=["studyId","sampleId","patientId"], columns="geneSymbol", values="value").reset_index()

merged = expr_wide.merge(surv, on=["studyId","patientId"], how="inner")
kirc = merged[merged["studyId"] == "kirc_tcga_pan_can_atlas_2018"][["time","event","RELA"]].dropna()

print(f"KIRC/RELA analysis set: n={len(kirc)}, events={kirc['event'].sum()}")

# --- Median split ---
median_val = kirc["RELA"].median()
kirc["group"] = np.where(kirc["RELA"] >= median_val, "High RELA", "Low RELA")
high = kirc[kirc["group"] == "High RELA"]
low = kirc[kirc["group"] == "Low RELA"]

# --- Standard, correct statistical test for comparing two KM curves: log-rank test ---
lr = logrank_test(high["time"], low["time"], high["event"], low["event"])

# --- Cox HR (continuous, matches the confirmed FDR-significant result) for effect size + CI ---
cph = CoxPHFitter()
cph.fit(kirc.rename(columns={"RELA": "expr"})[["time","event","expr"]], duration_col="time", event_col="event")
hr = cph.hazard_ratios_["expr"]
ci = cph.confidence_intervals_.loc["expr"]
hr_lo, hr_hi = np.exp(ci.iloc[0]), np.exp(ci.iloc[1])

print(f"Median RELA (z-score): {median_val:.3f}")
print(f"Log-rank test: p = {lr.p_value:.2e}")
print(f"Cox HR (continuous): {hr:.3f} (95% CI: {hr_lo:.3f}-{hr_hi:.3f})")

# --- Plot ---
fig, ax = plt.subplots(figsize=(7, 6))

kmf_high = KaplanMeierFitter()
kmf_high.fit(high["time"], high["event"], label=f"High RELA (n={len(high)})")
kmf_high.plot_survival_function(ax=ax, color="#C00000", ci_show=True)

kmf_low = KaplanMeierFitter()
kmf_low.fit(low["time"], low["event"], label=f"Low RELA (n={len(low)})")
kmf_low.plot_survival_function(ax=ax, color="#4472C4", ci_show=True)

ax.set_xlabel("Time (months)")
ax.set_ylabel("Overall survival probability")
ax.set_title("KIRC: Overall survival by RELA expression", fontsize=13)

# Annotate with the log-rank p-value and Cox HR -- the two numbers a reader needs
stat_text = f"Log-rank p = {lr.p_value:.2e}\nCox HR = {hr:.2f} (95% CI {hr_lo:.2f}-{hr_hi:.2f})"
ax.text(0.98, 0.05, stat_text, transform=ax.transAxes, ha="right", va="bottom",
        fontsize=10, bbox=dict(boxstyle="round", facecolor="white", edgecolor="gray"))

ax.legend(loc="upper right")
ax.set_ylim(0, 1.05)

# Standard number-at-risk table beneath the plot
add_at_risk_counts(kmf_high, kmf_low, ax=ax)

plt.tight_layout()
plt.savefig("KM_KIRC_RELA.png", dpi=300, bbox_inches="tight")
plt.savefig("KM_KIRC_RELA.eps", format="eps", bbox_inches="tight")
print("\nSaved KM_KIRC_RELA.png and .eps")

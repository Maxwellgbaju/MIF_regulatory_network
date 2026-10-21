"""Figure 4: pipeline dependence of RELA. Inputs: results/rev1b_rela_crosscheck.csv, results/rev1c_survival_by_pipeline.csv,
results/rev1d_agreement.csv, results/rev1d_nfkb_targets.csv."""
import numpy as np, pandas as pd
from figstyle import plt, save, tag, need
import matplotlib.ticker as mt
from matplotlib.patches import Patch

b1 = pd.read_csv(need("results/rev1b_rela_crosscheck.csv"))
c1 = pd.read_csv(need("results/rev1c_survival_by_pipeline.csv"))
ag = pd.read_csv(need("results/rev1d_agreement.csv"))
nf = pd.read_csv(need("results/rev1d_nfkb_targets.csv"))
SRC = ["Legacy_RSEM", "PanCanAtlas_EB++", "Toil_TPM"]
COL = {"Legacy_RSEM": "#7F7F7F", "Original_RSEM": "#7F7F7F", "PanCanAtlas_EB++": "#E69F00", "Toil_TPM": "#009E73"}
NAME = {"Legacy_RSEM": "Original TCGA RSEM", "PanCanAtlas_EB++": "PanCanAtlas (batch-corrected)", "Toil_TPM": "Toil re-quantification"}
cohorts = ["BRCA", "KIRC", "LIHC", "LUSC", "SKCM", "STAD"]

fig, axs = plt.subplots(2, 2, figsize=(7.2, 6.2), gridspec_kw=dict(hspace=0.6, wspace=0.45))
a = axs[0, 0]; xg, ticks = 0, []
for g in ["RELA", "SP1", "HIF1A"]:
    for i, c in enumerate(cohorts):
        for k, s in enumerate(SRC):
            v = b1[(b1.gene == g) & (b1.cohort == c) & (b1.source == s)].rho
            if len(v): a.scatter(xg + i + (k - 1) * 0.22, v.iloc[0], color=COL[s], s=11, zorder=3)
    ticks.append(xg + 2.5); xg += 7.5
a.axhline(0, color="k", lw=0.6); a.set_xticks(ticks); a.set_xticklabels(["RELA", "SP1", "HIF1A"])
a.set_ylabel("Spearman rho with MIF"); a.set_title("Co-expression with MIF, six cohorts", fontsize=8.3); tag(a, "a")
a.legend(handles=[Patch(color=COL[s], label=NAME[s]) for s in SRC], fontsize=5.8, frameon=False, loc="lower left")

b = axs[0, 1]; k_c = c1[c1.cohort == "KIRC"].copy()
lohi = k_c.CI.str.split("-", expand=True).astype(float); k_c["lo"], k_c["hi"] = lohi[0], lohi[1]
genes = ["RELA", "MIF", "SP1", "HIF1A"]
for gi, g in enumerate(genes):
    for k, s in enumerate(SRC):
        row = k_c[(k_c.gene == g) & (k_c.source == s)]
        if len(row):
            yy = gi * 4 + k; r = row.iloc[0]
            b.plot([r.lo, r.hi], [yy, yy], color=COL[s], lw=1.4); b.scatter(r.HR_per_SD, yy, color=COL[s], s=14, zorder=3)
b.axvline(1, color="k", lw=0.6, ls="--"); b.set_yticks([1, 5, 9, 13]); b.set_yticklabels(genes); b.invert_yaxis()
b.set_xscale("log"); b.set_xticks([0.7, 1, 1.4]); b.set_xticklabels(["0.7", "1.0", "1.4"]); b.xaxis.set_minor_formatter(mt.NullFormatter())
n_kirc = int(k_c.n.iloc[0]); ev = int(k_c.events.iloc[0])
b.set_xlabel("Overall survival HR per SD (95% CI)"); b.set_title(f"Kidney clear cell survival (n={n_kirc}, {ev} events)", fontsize=8.3); tag(b, "b")

c = axs[1, 0]; order = ["MIF", "HIF1A", "SP1", "RELA"]; coh3 = ["KIRC", "BRCA", "SKCM"]; w = 0.25
for k, (coh, colr) in enumerate(zip(coh3, ["#4C72B0", "#55A868", "#C44E52"])):
    vals = [ag[(ag.cohort == coh) & (ag.gene == g)].rho_toil_vs_original.iloc[0] for g in order]
    c.bar(np.arange(4) + (k - 1) * w, vals, width=w, color=colr, label=coh)
c.set_xticks(range(4)); c.set_xticklabels(order); c.set_ylim(0, 1.18); c.set_yticks([0, .2, .4, .6, .8, 1.0])
c.set_ylabel("Patient-level agreement\n(Spearman, Toil vs original)"); c.legend(fontsize=6.3, frameon=False, ncol=3, loc="upper center")
c.set_title("Same patients, two pipelines", fontsize=8.3); tag(c, "c")

d = axs[1, 1]; x = np.arange(3)
med = nf.groupby(["cohort", "source"]).rho_RELA_vs_target.median()
nk = nf[nf.target == "NFKB1"].set_index(["cohort", "source"]).rho_RELA_vs_target
o = [med[(k, "Original_RSEM")] for k in coh3]; t_ = [med[(k, "Toil_TPM")] for k in coh3]
d.bar(x - 0.2, o, 0.38, color=COL["Original_RSEM"], label="Original TCGA RSEM"); d.bar(x + 0.2, t_, 0.38, color=COL["Toil_TPM"], label="Toil")
for xi, k in zip(x, coh3):
    d.text(xi, max(med[(k, "Original_RSEM")], med[(k, "Toil_TPM")]) + 0.03,
           f"NFKB1: {nk[(k, 'Original_RSEM')]:.2f} vs {nk[(k, 'Toil_TPM')]:.2f}", ha="center", fontsize=6)
d.set_xticks(x); d.set_xticklabels(coh3); d.set_ylim(0, max(max(o), max(t_)) + 0.18)
d.set_ylabel("Median rho, RELA vs 7 NF-kB targets"); d.legend(fontsize=6.3, frameon=False, loc="upper right")
d.set_title("RELA vs its NF-kB target genes", fontsize=8.3); tag(d, "d")
save(fig, "Figure4_pipeline_dependence")

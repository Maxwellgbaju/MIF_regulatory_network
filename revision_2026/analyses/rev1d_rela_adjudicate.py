"""Updated rev1d: same analysis as before, now also SAVES results for Figure 4c-d.
Outputs: results/rev1d_agreement.csv, results/rev1d_nfkb_targets.csv"""
from mifrev_utils import *
from scipy.stats import spearmanr

CORE = ["MIF", "RELA", "SP1", "HIF1A"]
NFKB_TARGETS = ["NFKBIA", "TNFAIP3", "NFKB2", "RELB", "ICAM1", "NFKB1", "BIRC3"]
GENES = CORE + NFKB_TARGETS
COHORTS = {"Kidney Clear Cell Carcinoma": "KIRC", "Breast Invasive Carcinoma": "BRCA", "Skin Cutaneous Melanoma": "SKCM"}
meta = load_meta(); groups = tumor_groups(meta, 30)
all_t = sorted({s for v in groups.values() for s in v})
ids, miss = resolve(GENES)
toil = fetch(ids, all_t, "adjudicate"); toil.index = [s[:15] for s in toil.index]
toil = toil[~toil.index.duplicated()]
LEG, PAN = "https://tcga.xenahubs.net", "https://pancanatlas.xenahubs.net"
EB = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
pss = xena.dataset_samples(PAN, EB, None)
v = xena.dataset_probe_values(PAN, EB, pss, GENES)
pan = pd.DataFrame({g: pd.to_numeric(pd.Series(x, index=pss), errors="coerce") for g, x in zip(GENES, v[1])})
pan.index = [s[:15] for s in pan.index]; pan = pan[~pan.index.duplicated()]

agree, nfkb = [], []
for dis, code in COHORTS.items():
    ds = f"TCGA.{code}.sampleMap/HiSeqV2"
    ss = [s for s in xena.dataset_samples(LEG, ds, None) if s[13:15] == "01"]
    v = xena.dataset_probe_values(LEG, ds, ss, GENES)
    leg = pd.DataFrame({g: pd.to_numeric(pd.Series(x, index=ss), errors="coerce") for g, x in zip(GENES, v[1])})
    leg.index = [s[:15] for s in leg.index]; leg = leg[~leg.index.duplicated()]
    bar = sorted(set(s[:15] for s in groups[dis]) & set(leg.index))
    for g in CORE:
        d = pd.concat([toil.reindex(bar)[g], leg.reindex(bar)[g]], axis=1).dropna()
        agree.append({"cohort": code, "gene": g, "n": len(d), "rho_toil_vs_original": spearmanr(d.iloc[:, 0], d.iloc[:, 1])[0]})
    for label, e in [("Original_RSEM", leg), ("PanCanAtlas_EB++", pan), ("Toil_TPM", toil)]:
        sub = e.reindex(bar)
        for t in NFKB_TARGETS:
            if t in sub:
                nfkb.append({"cohort": code, "source": label, "target": t,
                             "rho_RELA_vs_target": spearmanr(sub["RELA"], sub[t], nan_policy="omit")[0]})
pd.DataFrame(agree).to_csv("results/rev1d_agreement.csv", index=False)
pd.DataFrame(nfkb).to_csv("results/rev1d_nfkb_targets.csv", index=False)
print(pd.DataFrame(agree).pivot(index="gene", columns="cohort", values="rho_toil_vs_original").round(3))
print(pd.DataFrame(nfkb).groupby(["cohort", "source"]).rho_RELA_vs_target.median().round(3))
print("saved results/rev1d_agreement.csv and results/rev1d_nfkb_targets.csv")

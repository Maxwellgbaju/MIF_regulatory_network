from mifrev_utils import *
from lifelines import CoxPHFitter

GENES = ["MIF", "RELA", "SP1", "HIF1A"]
COHORTS = {"Kidney Clear Cell Carcinoma": "KIRC", "Breast Invasive Carcinoma": "BRCA",
           "Liver Hepatocellular Carcinoma": "LIHC", "Stomach Adenocarcinoma": "STAD",
           "Lung Squamous Cell Carcinoma": "LUSC", "Skin Cutaneous Melanoma": "SKCM"}
meta = load_meta()
groups = tumor_groups(meta, 30)
all_t = sorted({s for v in groups.values() for s in v})

# survival (same dataset as the September Xena analysis), indexed by 15-char barcode
sv = xena.dataset_probe_values(TOIL, "TCGA_survival_data", all_t, ["OS", "OS.time"])
surv = pd.DataFrame({"OS": pd.to_numeric(pd.Series(sv[1][0], index=all_t), errors="coerce"),
                     "OS_time": pd.to_numeric(pd.Series(sv[1][1], index=all_t), errors="coerce")})
surv.index = [s[:15] for s in surv.index]
surv = surv[~surv.index.duplicated()]

# three expression sources
ids, _ = resolve(GENES)
toil = fetch(ids, all_t, "crosscheck"); toil.index = [s[:15] for s in toil.index]

LEG = "https://tcga.xenahubs.net"
def legacy(code):
    ds = f"TCGA.{code}.sampleMap/HiSeqV2"
    ss = [s for s in xena.dataset_samples(LEG, ds, None) if s[13:15] == "01"]
    v = xena.dataset_probe_values(LEG, ds, ss, GENES)
    return pd.DataFrame({g: pd.to_numeric(pd.Series(x, index=ss), errors="coerce") for g, x in zip(GENES, v[1])})

PAN = "https://pancanatlas.xenahubs.net"
EB = "EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena"
pss = xena.dataset_samples(PAN, EB, None)
v = xena.dataset_probe_values(PAN, EB, pss, GENES)
pan = pd.DataFrame({g: pd.to_numeric(pd.Series(x, index=pss), errors="coerce") for g, x in zip(GENES, v[1])})
pan.index = [s[:15] for s in pan.index]
pan = pan[~pan.index.duplicated()]

rows = []
for dis, code in COHORTS.items():
    bar = sorted({s[:15] for s in groups[dis]})
    sources = {"Toil_TPM": toil, "PanCanAtlas_EB++": pan}
    try:
        lg = legacy(code); lg.index = [s[:15] for s in lg.index]
        sources["Legacy_RSEM"] = lg[~lg.index.duplicated()]
    except Exception as ex:
        print(f"legacy {code} failed: {ex}")
    for src, e in sources.items():
        for g in GENES:
            if g not in e:
                continue
            d = pd.concat([e.reindex(bar)[g], surv.reindex(bar)], axis=1).dropna()
            d = d[d.OS_time > 0]
            if len(d) < 50 or d.OS.sum() < 10:
                continue
            d["z"] = (d[g] - d[g].mean()) / d[g].std()
            cph = CoxPHFitter().fit(d[["z", "OS_time", "OS"]], "OS_time", "OS")
            s = cph.summary.loc["z"]
            rows.append({"cohort": code, "gene": g, "source": src, "n": len(d), "events": int(d.OS.sum()),
                         "HR_per_SD": round(s["exp(coef)"], 3),
                         "CI": f"{s['exp(coef) lower 95%']:.2f}-{s['exp(coef) upper 95%']:.2f}",
                         "p": s["p"]})
res = pd.DataFrame(rows)
res["FDR_within_source"] = res.groupby("source")["p"].transform(lambda x: bh(x))
res.to_csv("results/rev1c_survival_by_pipeline.csv", index=False)
print("=== Overall survival HR per SD of expression, same patients, three pipelines ===")
for g in GENES:
    print(f"\n--- {g} ---")
    t = res[res.gene == g].pivot_table(index="cohort", columns="source", values="HR_per_SD")
    print(t.round(3).to_string())
print("\nKIRC detail:")
print(res[res.cohort == "KIRC"].sort_values(["gene", "source"]).to_string(index=False))
print("\nManuscript reference: RELA in KIRC HR = 1.31 (95% CI 1.14-1.51), FDR = 0.0035")

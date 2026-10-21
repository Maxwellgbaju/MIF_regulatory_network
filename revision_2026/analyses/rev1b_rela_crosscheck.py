from mifrev_utils import *
from scipy.stats import spearmanr

GENES = ["MIF", "RELA", "SP1", "HIF1A"]
COHORTS = {"Kidney Clear Cell Carcinoma": "KIRC", "Breast Invasive Carcinoma": "BRCA",
           "Liver Hepatocellular Carcinoma": "LIHC", "Stomach Adenocarcinoma": "STAD",
           "Lung Squamous Cell Carcinoma": "LUSC", "Skin Cutaneous Melanoma": "SKCM"}
meta = load_meta()
groups = tumor_groups(meta, 30)

def corr_table(get_expr, label):
    out = []
    for dis, code in COHORTS.items():
        e = get_expr(dis, code)
        if e is None or len(e) < 30:
            continue
        for g in ["RELA", "SP1", "HIF1A"]:
            d = e[["MIF", g]].dropna()
            rho, p = spearmanr(d.MIF, d[g])
            out.append({"source": label, "cohort": code, "gene": g, "n": len(d), "rho": round(rho, 3), "p": p})
    return out

rows = []
# 1. Toil (already cached)
ids, _ = resolve(GENES)
toil = fetch(ids, sorted({s for v in groups.values() for s in v}), "crosscheck")
rows += corr_table(lambda dis, code: toil.loc[groups[dis]], "Toil_TPM")

# 2. Original TCGA RSEM (legacy hub, per-cohort HiSeqV2, gene symbols)
LEG = "https://tcga.xenahubs.net"
def legacy(dis, code):
    ds = f"TCGA.{code}.sampleMap/HiSeqV2"
    try:
        ss = [s for s in xena.dataset_samples(LEG, ds, None) if s[13:15] == "01"]
        v = xena.dataset_probe_values(LEG, ds, ss, GENES)
        return pd.DataFrame({g: pd.to_numeric(pd.Series(x, index=ss), errors="coerce") for g, x in zip(GENES, v[1])})
    except Exception as ex:
        print(f"legacy {code} failed: {ex}"); return None
rows += corr_table(legacy, "Legacy_RSEM")

# 3. PanCanAtlas batch-corrected matrix (basis of cBioPortal PanCancer Atlas studies)
PAN = "https://pancanatlas.xenahubs.net"
try:
    dsets = xena.dataset_list(PAN, xena.all_cohorts(PAN, []))
except Exception:
    dsets = xena.dataset_list(PAN, [])
eb = next((d["name"] for d in dsets if "EB++" in d["name"] and "geneExp" in d["name"]), None)
print("PanCanAtlas expression dataset:", eb)
if eb:
    pss = xena.dataset_samples(PAN, eb, None)
    v = xena.dataset_probe_values(PAN, eb, pss, GENES)
    pan = pd.DataFrame({g: pd.to_numeric(pd.Series(x, index=pss), errors="coerce") for g, x in zip(GENES, v[1])})
    pan.index = [s[:15] for s in pan.index]
    def pancan(dis, code):
        ss = [s[:15] for s in groups[dis] if s[:15] in pan.index]
        return pan.loc[ss]
    rows += corr_table(pancan, "PanCanAtlas_EB++")

res = pd.DataFrame(rows)
res.to_csv("results/rev1b_rela_crosscheck.csv", index=False)
print("\n=== Spearman rho with MIF, same cohorts, three processing pipelines ===")
print(res.pivot_table(index=["gene", "cohort"], columns="source", values="rho").round(3).to_string())
print("\nSeptember cBioPortal values for reference: RELA KIRC +0.316, BRCA +0.229, SKCM +0.173, STAD +0.129, LIHC +0.078, LUSC -0.012")

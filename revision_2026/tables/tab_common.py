"""Shared helpers for the manuscript tables. Run from /shared/rc/malariarna/MIF_revision_2026 (needs mifrev_utils.py)."""
import os
import numpy as np, pandas as pd
from scipy.stats import spearmanr, pearsonr
from mifrev_utils import xena, load_meta, tumor_groups, resolve, fetch, bh

os.makedirs("tables", exist_ok=True)
LEG = "https://tcga.xenahubs.net"
COHORTS = {"Kidney Clear Cell Carcinoma": "KIRC", "Lung Adenocarcinoma": "LUAD", "Liver Hepatocellular Carcinoma": "LIHC",
           "Stomach Adenocarcinoma": "STAD", "Breast Invasive Carcinoma": "BRCA", "Lung Squamous Cell Carcinoma": "LUSC",
           "Skin Cutaneous Melanoma": "SKCM"}
STROMA = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRB"]
IMMUNE = ["PTPRC", "CD3E", "CD2", "LCP2", "CD53"]

def need(path):
    if not os.path.exists(path):
        raise SystemExit(f"Missing input: {path}")
    return path

def legacy_expr(code, genes):
    """Original TCGA RSEM quantification (Xena TCGA hub), primary tumours, indexed by 15-character barcode."""
    ds = f"TCGA.{code}.sampleMap/HiSeqV2"
    ss = [s for s in xena.dataset_samples(LEG, ds, None) if s[13:15] == "01"]
    v = xena.dataset_probe_values(LEG, ds, ss, genes)
    df = pd.DataFrame({g: pd.to_numeric(pd.Series(x, index=ss), errors="coerce") for g, x in zip(genes, v[1])})
    df.index = [s[:15] for s in df.index]
    return df[~df.index.duplicated()]

def toil_expr(genes, tag):
    meta = load_meta(); groups = tumor_groups(meta, 30)
    samples = sorted({s for v in groups.values() for s in v})
    ids, missing = resolve(genes)
    if missing: print("  unresolved in Toil:", missing)
    df = fetch(ids, samples, tag)
    df.index = [s[:15] for s in df.index]
    df = df[~df.index.duplicated()]
    bars = {dis: sorted({s[:15] for s in groups.get(dis, [])}) for dis in COHORTS}
    return df, bars

def zscore(df):
    sd = df.std(ddof=0)
    return ((df - df.mean()) / sd).loc[:, sd > 0]

def composition_score(e):
    z = zscore(e)
    st = z[[g for g in STROMA if g in z]].mean(axis=1)
    im = z[[g for g in IMMUNE if g in z]].mean(axis=1)
    return st, im

def partial_spearman(x, y, covs):
    rx, ry = pd.Series(x).rank().values, pd.Series(y).rank().values
    C = np.column_stack([pd.Series(c).rank().values for c in covs] + [np.ones(len(rx))])
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    r = np.corrcoef(ex, ey)[0, 1]
    df = len(rx) - 2 - len(covs)
    from scipy.stats import t as tdist
    tt = r * np.sqrt(df / max(1e-12, 1 - r ** 2))
    return r, 2 * tdist.sf(abs(tt), df)

def stars(fdr):
    return "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*" if fdr < 0.05 else ""

def show(title, df, index=True):
    """Print a finished table as plain CSV between markers, for copying."""
    print(f"\n===== BEGIN {title} (copy from here) =====")
    print(df.to_csv(index=index).strip())
    print(f"===== END {title} =====\n")

"""Table 1: integrated evidence for 16 candidate transcription factors.
Inputs: results/rev1_MIF_coexpression_all.csv, results/rev1b_rela_crosscheck.csv, results/rev2_perturbation_all_terms.csv,
results/rev3_tf_activity_vs_MIF.csv; ChIP-Atlas target tables at 1, 5 and 10 kb (downloaded/cached here).
Output: tables/Table1_integrated_evidence.csv"""
import os, requests
import numpy as np, pandas as pd
from tab_common import need, show

TFS = ["HIF1A", "TFAP2A", "SP1", "RELA", "NFKB1", "MYC", "TP53", "ESR1", "EGR1", "KLF4",
       "IRF1", "IRF2", "NFIL3", "NR1H2", "RXRA", "REST"]
# Nomination evidence and maximum JASPAR score, as reported in the manuscript (evidence tiers; Figure 3 motif screen)
NOMINATION = {
    "SP1": ("Motif, conserved, TRRUST, ENCODE (4 of 4 tiers)", 17.55), "HIF1A": ("Motif, conserved, TRRUST (3 of 4)", None),
    "RELA": ("Motif, conserved, TRRUST (3 of 4)", None), "NFKB1": ("TRRUST; below conservation threshold", None),
    "EGR1": ("Hypothesis-free motif screen", 18.20), "NR1H2": ("Hypothesis-free screen (NR1H2::RXRA)", 16.51),
    "RXRA": ("Hypothesis-free screen (NR1H2::RXRA)", 16.51), "IRF1": ("Hypothesis-free motif screen", 15.74),
    "KLF4": ("Hypothesis-free motif screen", 15.00), "IRF2": ("Hypothesis-free motif screen", 14.97),
    "NFIL3": ("Hypothesis-free motif screen", 14.92), "ESR1": ("Hypothesis-free motif screen", 14.79),
    "TP53": ("Hypothesis-free motif screen", 14.70), "TFAP2A": ("Hypothesis-free motif screen", 14.57),
    "MYC": ("Hypothesis-free motif screen", 14.48), "REST": ("Hypothesis-free motif screen", 14.36)}

def chipatlas(tf, dist):
    p = f"cache/chipatlas_{tf}_1kb.tsv" if dist == 1 else f"cache/chipatlas_{tf}_{dist}kb.tsv"
    if not os.path.exists(p):
        r = requests.get(f"https://chip-atlas.dbcls.jp/data/hg38/target/{tf}.{dist}.tsv", timeout=1800)
        if r.status_code != 200:
            return None
        open(p, "wb").write(r.content)
    df = pd.read_csv(p, sep="\t")
    gene_col, avg_col = df.columns[0], df.columns[1]
    exp = [c for c in df.columns[2:] if "STRING" not in c.upper()]
    row = df[df[gene_col] == "MIF"]
    if row.empty:
        return {"n": len(exp), "bound": 0, "pct": np.nan}
    v = row[exp].iloc[0].astype(float)
    return {"n": len(exp), "bound": int((v > 0).sum()), "pct": 100 * (df[avg_col] < row[avg_col].iloc[0]).mean()}

co = pd.read_csv(need("results/rev1_MIF_coexpression_all.csv"))
cb = pd.read_csv(need("results/rev1b_rela_crosscheck.csv"))
pt = pd.read_csv(need("results/rev2_perturbation_all_terms.csv"))
ac = pd.read_csv(need("results/rev3_tf_activity_vs_MIF.csv"))

rows = []
for tf in TFS:
    print(f"{tf} ...", flush=True)
    nom, score = NOMINATION[tf]
    r = {"Factor": tf, "Nomination evidence": nom, "Max JASPAR score": "" if score is None else f"{score:.2f}"}
    for d in (1, 5, 10):
        c = chipatlas(tf, d)
        r[f"ChIP-Atlas {d} kb: bound/total"] = "n.a." if c is None else f"{c['bound']}/{c['n']}"
        r[f"ChIP-Atlas {d} kb: percentile"] = "n.a." if c is None or np.isnan(c["pct"]) else f"{c['pct']:.1f}"
    p = pt[pt.TF == tf]
    act = (p.interpretation == "consistent with TF activating MIF").sum()
    rep = (p.interpretation == "consistent with TF repressing MIF").sum()
    r["Perturbation signatures (MIF responded)"] = f"{len(p)} ({int(p.MIF_in_set.sum())})"
    r["Directional: activating / repressing"] = f"{act} / {rep}"
    c1 = co[co.gene == tf]
    r["Co-expression, +/- cancer types (median rho)"] = ("not tested" if c1.empty else
        f"{int(((c1.rho > 0) & (c1.FDR < 0.05)).sum())} / {int(((c1.rho < 0) & (c1.FDR < 0.05)).sum())} ({c1.rho.median():+.2f})")
    a1 = ac[ac.TF == tf]
    r["Activity vs MIF, median rho (significant)"] = ("not tested" if a1.empty else
        f"{a1.rho_activity_vs_MIF.median():+.2f} ({int((a1.FDR_activity_perm < 0.05).sum())})")
    b1 = cb[cb.gene == tf]
    if b1.empty:
        r["Direction agrees across pipelines (6 cohorts)"] = "not tested"
    else:
        w = b1.pivot_table(index="cohort", columns="source", values="rho")
        agree = int((np.sign(w["Toil_TPM"]) == np.sign(w["Legacy_RSEM"])).sum())
        r["Direction agrees across pipelines (6 cohorts)"] = f"{agree} of {len(w)}"
    rows.append(r)

out = pd.DataFrame(rows)
out.to_csv("tables/Table1_integrated_evidence.csv", index=False)
show("TABLE 1", out, index=False)

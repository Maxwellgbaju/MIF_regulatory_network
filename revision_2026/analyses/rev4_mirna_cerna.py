import sys
from mifrev_utils import *
from scipy.stats import spearmanr

PAN = "https://pancanatlas.xenahubs.net"
MIRNAS = ["hsa-miR-92a-1-5p", "hsa-miR-654-5p", "hsa-miR-608", "hsa-miR-541-3p", "hsa-miR-451a",
          "hsa-miR-363-5p", "hsa-miR-23b-5p", "hsa-miR-23a-5p", "hsa-miR-1537-5p", "hsa-miR-1228-5p",
          "hsa-miR-744-5p", "hsa-miR-320a", "hsa-miR-146a-5p", "hsa-miR-1246", "hsa-miR-769-5p"]
LNCRNAS = ["AL356488.2", "NORAD", "LMNTD2-AS1", "AL355987.4", "AC020922.3", "MALAT1", "MIR3142HG"]
PAIRS = [(m, l) for m in ["hsa-miR-654-5p", "hsa-miR-541-3p"]
         for l in ["AL356488.2", "LMNTD2-AS1", "AL355987.4", "AC020922.3"]] + [
        ("hsa-miR-744-5p", "AL356488.2"), ("hsa-miR-451a", "NORAD"), ("hsa-miR-146a-5p", "MALAT1"),
        ("hsa-miR-769-5p", "MALAT1"), ("hsa-miR-1246", "MALAT1"), ("hsa-miR-320a", "MIR3142HG")]

try:
    dsets = xena.dataset_list(PAN, xena.all_cohorts(PAN, []))
except Exception:
    dsets = xena.dataset_list(PAN, [])
names = [d["name"] for d in dsets if "mir" in d["name"].lower()]
print("miRNA datasets on PanCanAtlas hub:", names)
ds = sys.argv[1] if len(sys.argv) > 1 else next((n for n in names if "pancanmirs" in n.lower()), names[0] if names else None)
if ds is None:
    raise SystemExit("No miRNA dataset found; share the list above.")
print("Using:", ds)

fields = xena.dataset_field(PAN, ds)
low = {f.lower(): f for f in fields}
matched = {}
for m in MIRNAS:
    if m.lower() in low:
        matched[m] = low[m.lower()]
    else:
        core = m.lower().replace("hsa-", "")
        cand = [f for f in fields if core in f.lower()]
        if cand:
            matched[m] = cand[0]
print("\nmiRNA name matches (check these):")
for m in MIRNAS:
    print(f"  {m:20s} -> {matched.get(m, 'NOT FOUND')}")

msamples = xena.dataset_samples(PAN, ds, None)
mres = xena.dataset_probe_values(PAN, ds, msamples, list(matched.values()))
mir = pd.DataFrame({m: pd.to_numeric(pd.Series(v, index=msamples), errors="coerce")
                    for m, v in zip(matched, mres[1])})
mir.index = [s[:15] for s in mir.index]
mir = mir[~mir.index.duplicated()]

meta = load_meta()
groups = tumor_groups(meta, 30)
samples = sorted({s for v in groups.values() for s in v})
ids, missing = resolve(["MIF"] + LNCRNAS)
print("\nUnresolved lncRNAs (clone-based names often renamed):", missing)
expr = fetch(ids, samples, "cerna")
expr.index = [s[:15] for s in expr.index]
expr = expr[~expr.index.duplicated()]

rows = []
for dis, ss in groups.items():
    ss = sorted({s[:15] for s in ss if s[:15] in mir.index})
    if len(ss) < 30:
        continue
    e, mm = expr.reindex(ss), mir.loc[ss]
    def test(kind, a, b, x, y, expected):
        d = pd.concat([x, y], axis=1).dropna()
        if len(d) < 30 or d.iloc[:, 0].nunique() < 10 or d.iloc[:, 1].nunique() < 10:
            return
        rho, p = spearmanr(d.iloc[:, 0], d.iloc[:, 1])
        rows.append({"test": kind, "a": a, "b": b, "cancer": dis, "n": len(d), "rho": rho,
                     "p": p, "expected": expected})
    for m in matched:
        test("miRNA_vs_MIF", m, "MIF", mm[m], e["MIF"], "negative")
    for l in ids:
        if l != "MIF":
            test("lncRNA_vs_MIF", l, "MIF", e[l], e["MIF"], "positive")
    for m, l in PAIRS:
        if m in matched and l in ids:
            test("miRNA_vs_lncRNA", m, l, mm[m], e[l], "negative")

res = pd.DataFrame(rows)
res["FDR"] = bh(res.p)
res["sig_expected"] = (res.FDR < 0.05) & (((res.expected == "negative") & (res.rho < 0)) | ((res.expected == "positive") & (res.rho > 0)))
res["sig_opposite"] = (res.FDR < 0.05) & ~res.sig_expected
res.to_csv("results/rev4_mirna_cerna_all.csv", index=False)
summ = res.groupby(["test", "a", "b"]).agg(cancers=("cancer", "count"), median_rho=("rho", "median"),
        sig_expected_dir=("sig_expected", "sum"), sig_opposite_dir=("sig_opposite", "sum")).reset_index()
summ.to_csv("results/rev4_mirna_cerna_summary.csv", index=False)
print("\n=== Expression test of the sponge network (expected: miRNA-MIF negative, lncRNA-MIF positive, miRNA-lncRNA negative) ===")
print(summ.round(3).to_string(index=False))

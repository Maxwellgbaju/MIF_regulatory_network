from mifrev_utils import *
import gseapy as gp
from scipy.stats import spearmanr

hall = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
PROGRAMS = ["Oxidative Phosphorylation", "Myc Targets V1", "Myc Targets V2", "Hypoxia", "Glycolysis"]
STROMA = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRB"]
IMMUNE = ["PTPRC", "CD3E", "CD2", "LCP2", "CD53"]
meta = load_meta(); groups = tumor_groups(meta, 30)
samples = sorted({s for v in groups.values() for s in v})

expr = pd.read_pickle("cache/expr_hallmark.pkl").loc[samples]
ids, _ = resolve(["MYC", "MAX"])
extra = fetch(ids, samples, "myc")
expr = pd.concat([expr, extra[[c for c in extra.columns if c not in expr.columns]]], axis=1)

prog_genes = {p: set(hall[p]) - {"MIF"} for p in PROGRAMS}
excluded = set().union(*prog_genes.values()) | set(STROMA) | set(IMMUNE) | {"MIF", "MYC", "MAX"}
means = expr.mean()
pool = [g for g in expr.columns if g not in excluded and expr[g].notna().mean() > 0.9]
controls = list((means[pool] - means["MIF"]).abs().sort_values().index[:100])
print(f"MIF mean log2 TPM {means['MIF']:.2f}; 100 expression-matched control genes, mean range "
      f"{means[controls].min():.2f}-{means[controls].max():.2f}", flush=True)

def partial(x, y, covs):
    rx, ry = pd.Series(x).rank().values, pd.Series(y).rank().values
    C = np.column_stack([pd.Series(c).rank().values for c in covs] + [np.ones(len(rx))])
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    return np.corrcoef(ex, ey)[0, 1]

rows, myc_rows = [], []
for dis, ss in groups.items():
    e = expr.loc[ss].dropna(axis=1, how="all"); sd = e.std(ddof=0)
    z = ((e - e.mean()) / sd).loc[:, sd > 0]
    covs = [z[[g for g in STROMA if g in z]].mean(axis=1), z[[g for g in IMMUNE if g in z]].mean(axis=1)]
    ok = e["MIF"].notna()
    covs = [c[ok] for c in covs]
    for p in PROGRAMS:
        sc = z[[g for g in prog_genes[p] if g in z]].mean(axis=1)[ok]
        mif_r = partial(sc, e["MIF"][ok], covs)
        ctrl_r = [partial(sc, e[c][ok].fillna(e[c].median()), covs) for c in controls if c in e]
        rows.append({"cancer": dis, "program": p, "MIF_partial_rho": mif_r,
                     "control_median_rho": np.median(ctrl_r),
                     "MIF_percentile_vs_controls": 100 * np.mean(np.array(ctrl_r) < mif_r)})
    if "MYC" in e:
        d = pd.concat([e["MYC"], e["MIF"]], axis=1).dropna()
        cv = [c.reindex(d.index) for c in [z[[g for g in STROMA if g in z]].mean(axis=1), z[[g for g in IMMUNE if g in z]].mean(axis=1)]]
        myc_rows.append({"cancer": dis, "n": len(d), "rho_MYC_MIF": spearmanr(d.MYC, d.MIF)[0],
                         "partial_rho_adj": partial(d.MYC, d.MIF, cv)})

res = pd.DataFrame(rows); res.to_csv("results/rev8_program_specificity.csv", index=False)
summ = res.groupby("program").agg(median_MIF_rho=("MIF_partial_rho", "median"),
                                  median_control_rho=("control_median_rho", "median"),
                                  median_MIF_percentile=("MIF_percentile_vs_controls", "median"),
                                  cancers_MIF_above_90th=("MIF_percentile_vs_controls", lambda x: int((x >= 90).sum())))
print("\n=== Is MIF special? Program association for MIF vs 100 expression-matched genes (composition-adjusted) ===")
print(summ.round(3).to_string())
m = pd.DataFrame(myc_rows); m.to_csv("results/rev8_MYC_vs_MIF.csv", index=False)
print("\n=== MYC mRNA vs MIF (Toil), composition-adjusted ===")
print(f"median raw rho {m.rho_MYC_MIF.median():.3f}; median adjusted rho {m.partial_rho_adj.median():.3f}; "
      f"adjusted positive in {(m.partial_rho_adj > 0).sum()} of {len(m)} cancers")
print(m.sort_values("partial_rho_adj", ascending=False).round(3).to_string(index=False))

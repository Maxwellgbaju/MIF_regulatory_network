from mifrev_utils import *
import gseapy as gp

hall = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
STROMA = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRB"]
IMMUNE = ["PTPRC", "CD3E", "CD2", "LCP2", "CD53"]
meta = load_meta(); groups = tumor_groups(meta, 30)
samples = sorted({s for v in groups.values() for s in v})
genes = sorted({g for v in hall.values() for g in v} | {"MIF"} | set(STROMA) | set(IMMUNE))
ids, missing = resolve(genes)
print(f"{len(hall)} Hallmark sets, {len(ids)} genes resolved ({len(missing)} unresolved)", flush=True)
expr = fetch(ids, samples, "hallmark")

def partial(x, y, covs):
    rx, ry = pd.Series(x).rank().values, pd.Series(y).rank().values
    C = np.column_stack([pd.Series(c).rank().values for c in covs] + [np.ones(len(rx))])
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    return np.corrcoef(ex, ey)[0, 1]

rows = []
for dis, ss in groups.items():
    e = expr.loc[ss]; sd = e.std(ddof=0)
    z = ((e - e.mean()) / sd).loc[:, sd > 0]
    mif = e["MIF"]
    stroma = z[[g for g in STROMA if g in z]].mean(axis=1)
    immune = z[[g for g in IMMUNE if g in z]].mean(axis=1)
    for name, gs in hall.items():
        g = [x for x in gs if x in z.columns and x != "MIF"]
        if len(g) < 15:
            continue
        sc = z[g].mean(axis=1)
        ok = ~(mif.isna() | sc.isna())
        rows.append({"cancer": dis, "set": name, "n_genes": len(g),
                     "partial_rho": partial(sc[ok], mif[ok], [stroma[ok], immune[ok]])})
res = pd.DataFrame(rows)
res["rank_in_cancer"] = res.groupby("cancer")["partial_rho"].rank(ascending=False)
res.to_csv("results/rev7_hallmark_vs_MIF_all.csv", index=False)
summ = res.groupby("set").agg(median_partial_rho=("partial_rho", "median"),
                              median_rank=("rank_in_cancer", "median"),
                              cancers_in_top5=("rank_in_cancer", lambda r: int((r <= 5).sum())),
                              cancers=("cancer", "count")).sort_values("median_rank")
summ.to_csv("results/rev7_hallmark_specificity_summary.csv")
print("\n=== Hallmark programs ranked by composition-adjusted association with MIF (rank 1 = strongest) ===")
print(summ.round(3).head(15).to_string())
print("\n... bottom 5:")
print(summ.round(3).tail(5).to_string())
print("\nHypoxia:", summ.loc["Hypoxia"].round(3).to_dict(), f"| overall position {list(summ.index).index('Hypoxia') + 1} of {len(summ)}")

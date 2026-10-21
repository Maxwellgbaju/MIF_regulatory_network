from mifrev_utils import *
import gseapy as gp
from scipy.stats import spearmanr

rng = np.random.default_rng(7)
NPERM = 200
hall = gp.get_library(name="MSigDB_Hallmark_2020", organism="Human")
HYP = [g for g in hall["Hypoxia"] if g != "MIF"]
STROMA = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FAP", "PDGFRB"]
IMMUNE = ["PTPRC", "CD3E", "CD2", "LCP2", "CD53"]

meta = load_meta(); groups = tumor_groups(meta, 30)
samples = sorted({s for v in groups.values() for s in v})
ids, missing = resolve(["MIF"] + HYP + STROMA + IMMUNE)
print(f"Hypoxia genes (MIF removed): {len(HYP)}; unresolved: {len(missing)}")
expr = fetch(ids, samples, "hypoxia")
pool_df = pd.read_pickle("cache/expr_tfactivity.pkl")          # background genes from Analysis 3
pool_df = pool_df[[c for c in pool_df.columns if c not in set(HYP) | {"MIF"} | set(STROMA) | set(IMMUNE)]]

def zmean(df):
    sd = df.std(ddof=0); z = (df - df.mean()) / sd
    return z.loc[:, sd > 0].mean(axis=1)

def partial(x, y, covs):
    rx, ry = pd.Series(x).rank().values, pd.Series(y).rank().values
    C = np.column_stack([pd.Series(c).rank().values for c in covs] + [np.ones(len(rx))])
    ex = rx - C @ np.linalg.lstsq(C, rx, rcond=None)[0]
    ey = ry - C @ np.linalg.lstsq(C, ry, rcond=None)[0]
    return np.corrcoef(ex, ey)[0, 1]

rows = []
for dis, ss in groups.items():
    e, pdf = expr.loc[ss], pool_df.loc[ss]
    mif = e["MIF"].values
    hyp = zmean(e[[g for g in HYP if g in e]]).values
    stroma = zmean(e[[g for g in STROMA if g in e]]).values
    immune = zmean(e[[g for g in IMMUNE if g in e]]).values
    ok = ~np.isnan(mif) & ~np.isnan(hyp)
    rho = spearmanr(hyp[ok], mif[ok])[0]
    prho = partial(hyp[ok], mif[ok], [stroma[ok], immune[ok]])
    k = len([g for g in HYP if g in e])
    null_r, null_p = [], []
    for _ in range(NPERM):
        rs = rng.choice(pdf.columns, size=min(k, pdf.shape[1]), replace=False)
        sc = zmean(pdf[rs]).values
        null_r.append(spearmanr(sc[ok], mif[ok])[0])
        null_p.append(partial(sc[ok], mif[ok], [stroma[ok], immune[ok]]))
    rows.append({"cancer": dis, "n": int(ok.sum()), "rho_hypoxia_MIF": rho, "partial_rho_adj_stroma_immune": prho,
                 "null_partial_median": np.median(null_p),
                 "perm_p_partial": (np.sum(np.array(null_p) >= prho) + 1) / (NPERM + 1)})
res = pd.DataFrame(rows).sort_values("partial_rho_adj_stroma_immune", ascending=False)
res["FDR"] = bh(res.perm_p_partial)
res.to_csv("results/rev6_hypoxia_vs_MIF.csv", index=False)
print(res.round(3).to_string(index=False))
pos = ((res.FDR < 0.05) & (res.partial_rho_adj_stroma_immune > res.null_partial_median)).sum()
print(f"\nCancers where the hypoxia score is positively associated with MIF beyond random gene sets (adjusted, FDR<0.05): {pos} of {len(res)}")
print(f"Median adjusted rho: hypoxia {res.partial_rho_adj_stroma_immune.median():.3f} vs random sets {res.null_partial_median.median():.3f}")

import io
from mifrev_utils import *
from scipy.stats import spearmanr

TFS = ["SP1", "HIF1A", "RELA", "TFAP2A", "EGR1", "KLF4", "IRF1", "IRF2", "NFIL3",
       "ESR1", "TP53", "NR1H2", "RXRA"]
CAP, NPERM = 250, 200
rng = np.random.default_rng(42)

def get_network():
    p = os.path.join(CACHE, "tf_network.tsv")
    if os.path.exists(p):
        return pd.read_csv(p, sep="\t")
    for name, url in [
        ("CollecTRI", "https://omnipathdb.org/interactions?datasets=collectri&genesymbols=yes&fields=curation_effort"),
        ("DoRothEA_ABC", "https://omnipathdb.org/interactions?datasets=dorothea&dorothea_levels=A,B,C&genesymbols=yes&fields=curation_effort")]:
        try:
            r = requests.get(url, timeout=600); r.raise_for_status()
            df = pd.read_csv(io.StringIO(r.text), sep="\t")
            if len(df) > 1000 and "source_genesymbol" in df:
                df["network"] = name
                df.to_csv(p, sep="\t", index=False)
                print(f"Using {name}: {len(df)} interactions")
                return df
        except Exception as e:
            print(f"{name} failed: {e}")
    raise SystemExit("Could not download a TF-target network from OmniPath.")

net = get_network()
print("Network:", net["network"].iloc[0])
if "consensus_inhibition" in net:
    net["sign"] = np.where((net.consensus_inhibition == 1) & (net.consensus_stimulation == 0), -1, 1)
else:
    net["sign"] = np.where(net.get("is_inhibition", 0) == 1, -1, 1)

targets, signs = {}, {}
for tf in TFS:
    sel = net[net.source_genesymbol.astype(str).str.split("_").apply(lambda t: tf in t)]
    if "MIF" in set(sel.target_genesymbol):
        print(f"  NOTE: MIF is itself a curated target of {tf} in this network")
    sel = sel[~sel.target_genesymbol.isin(["MIF", tf])]
    if "curation_effort" in sel:
        sel = sel.sort_values("curation_effort", ascending=False)
    sel = sel.drop_duplicates("target_genesymbol").head(CAP)
    targets[tf] = list(sel.target_genesymbol)
    signs[tf] = dict(zip(sel.target_genesymbol, sel.sign))
    print(f"{tf}: {len(targets[tf])} targets used")

meta = load_meta()
groups = tumor_groups(meta, 30)
samples = sorted({s for v in groups.values() for s in v})
genes = sorted(set(["MIF"] + TFS + [t for v in targets.values() for t in v]))
ids, missing = resolve(genes)
print(f"Fetching {len(ids)} genes x {len(samples)} tumors ({len(missing)} unresolved)...")
expr = fetch(ids, samples, "tfactivity")

rows = []
for dis, ss in groups.items():
    sub = expr.loc[ss]
    sd = sub.std(ddof=0)
    z = ((sub - sub.mean()) / sd).loc[:, sd > 0]
    mif = sub["MIF"]
    pool = [g for g in z.columns if g not in TFS and g != "MIF"]
    for tf in TFS:
        tg = [t for t in targets[tf] if t in z.columns]
        if len(tg) < 10 or tf not in sub:
            continue
        w = np.array([signs[tf][t] for t in tg])
        act = (z[tg].fillna(0).values * w).sum(1) / np.sqrt(len(tg))
        rho_act, p_act = spearmanr(act, mif, nan_policy="omit")
        rho_mrna, p_mrna = spearmanr(sub[tf], mif, nan_policy="omit")
        null = []
        for _ in range(NPERM):
            rs = rng.choice(pool, size=len(tg), replace=False)
            a = (z[rs].fillna(0).values * rng.permutation(w)).sum(1) / np.sqrt(len(tg))
            null.append(spearmanr(a, mif, nan_policy="omit")[0])
        perm_p = (np.sum(np.abs(null) >= abs(rho_act)) + 1) / (NPERM + 1)
        rows.append({"cancer": dis, "TF": tf, "n_targets": len(tg), "rho_mRNA_vs_MIF": rho_mrna,
                     "p_mRNA": p_mrna, "rho_activity_vs_MIF": rho_act, "p_activity": p_act,
                     "perm_p_activity": perm_p})
res = pd.DataFrame(rows)
res["FDR_mRNA"] = bh(res.p_mRNA)
res["FDR_activity_perm"] = bh(res.perm_p_activity)
res.to_csv("results/rev3_tf_activity_vs_MIF.csv", index=False)

summ = res.groupby("TF").agg(
    cancers=("cancer", "count"),
    median_rho_mRNA=("rho_mRNA_vs_MIF", "median"),
    median_rho_activity=("rho_activity_vs_MIF", "median"),
    activity_pos_sig=("FDR_activity_perm", lambda x: int(((x < 0.05) & (res.loc[x.index, "rho_activity_vs_MIF"] > 0)).sum())),
    activity_neg_sig=("FDR_activity_perm", lambda x: int(((x < 0.05) & (res.loc[x.index, "rho_activity_vs_MIF"] < 0)).sum())))
summ.to_csv("results/rev3_tf_activity_summary.csv")
print("\n=== mRNA vs inferred activity: correlation with MIF (activity significance = permutation, BH) ===")
print(summ.round(3).to_string())
print("\n=== SP1 by cancer type (the paradox) ===")
print(res[res.TF == "SP1"][["cancer", "rho_mRNA_vs_MIF", "rho_activity_vs_MIF", "FDR_activity_perm"]].round(3).to_string(index=False))
print("\n=== RELA by cancer type (target-based, independent of RELA mRNA quantification) ===")
print(res[res.TF == "RELA"][["cancer", "rho_mRNA_vs_MIF", "rho_activity_vs_MIF", "FDR_activity_perm"]].round(3).to_string(index=False))

import os, requests
import pandas as pd

TFS = ["SP1", "RELA", "NFKB1", "HIF1A", "TFAP2A", "EGR1", "KLF4", "IRF1", "IRF2", "NFIL3",
       "ESR1", "TP53", "NR1H2", "RXRA", "MYC", "REST"]
rows = []
for tf in TFS:
    p = f"cache/chipatlas_{tf}_1kb.tsv"
    if not os.path.exists(p):
        url = f"https://chip-atlas.dbcls.jp/data/hg38/target/{tf}.1.tsv"
        print(f"Downloading {tf}...", flush=True)
        r = requests.get(url, timeout=900)
        if r.status_code != 200:
            print(f"{tf}: not available in ChIP-Atlas (HTTP {r.status_code})")
            rows.append({"TF": tf, "experiments_total": 0, "note": "no ChIP-Atlas target file"})
            continue
        open(p, "wb").write(r.content)
    df = pd.read_csv(p, sep="\t")
    gene_col, avg_col = df.columns[0], df.columns[1]
    exp_cols = [c for c in df.columns[2:] if "STRING" not in c.upper()]
    row = df[df[gene_col] == "MIF"]
    if row.empty:
        rows.append({"TF": tf, "experiments_total": len(exp_cols), "experiments_binding_MIF": 0,
                     "note": "no peak within 1 kb of MIF TSS in any experiment"})
        print(f"{tf}: MIF not bound in {len(exp_cols)} experiments")
        continue
    vals = row[exp_cols].iloc[0].astype(float)
    bound = vals[vals > 0]
    cells = sorted({c.split("|")[1] if "|" in c else c for c in bound.index})
    pct = 100 * (df[avg_col] < row[avg_col].iloc[0]).mean()
    rows.append({"TF": tf, "experiments_total": len(exp_cols), "experiments_binding_MIF": len(bound),
                 "fraction_binding": round(len(bound) / max(len(exp_cols), 1), 3),
                 "MIF_score_percentile_among_targets": round(pct, 1),
                 "cell_types": "; ".join(cells[:25]) + (" ..." if len(cells) > 25 else "")})
    print(f"{tf}: MIF bound in {len(bound)}/{len(exp_cols)} experiments; score percentile {pct:.1f}", flush=True)
res = pd.DataFrame(rows)
res.to_csv("results/rev5_chipatlas_MIF_promoter.csv", index=False)
print("\n", res.drop(columns=["cell_types"], errors="ignore").to_string(index=False))

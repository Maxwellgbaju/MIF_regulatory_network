"""Shared helpers for the MIF revision analyses (reuses the September Xena setup)."""
import os, json, time
import requests
import numpy as np
import pandas as pd
import xenaPython as xena

TOIL = "https://toil.xenahubs.net"
TPM = "TcgaTargetGtex_rsem_gene_tpm"          # log2(TPM+0.001), same dataset as the manuscript
OLD = "/shared/rc/malariarna/MIF_student_network_v2/gepia_test"
CACHE = "cache"
os.makedirs(CACHE, exist_ok=True)
os.makedirs("results", exist_ok=True)

# IDs verified in September (GENCODE v23 names differ for some lncRNAs)
OVERRIDES = {"MIF": "ENSG00000240972", "NORAD": "ENSG00000260032",
             "AL356488.2": "ENSG00000270066"}

def load_all_ids():
    p = os.path.join(OLD, "xena_all_gene_ids.txt")
    if not os.path.exists(p):
        p = os.path.join(CACHE, "xena_all_gene_ids.txt")
        if not os.path.exists(p):
            open(p, "w").write("\n".join(xena.dataset_field(TOIL, TPM)))
    ids = [l.strip() for l in open(p) if l.strip()]
    return {i.split(".")[0]: i for i in ids}

def load_meta():
    p = os.path.join(OLD, "xena_phenotype_full.json")
    if not os.path.exists(p):
        raise SystemExit(f"Missing {p}: copy it from the September working folder.")
    return json.load(open(p))

def tumor_groups(meta, min_n=30):
    """TCGA primary tumors grouped by cancer type (same grouping as the manuscript)."""
    g = {}
    for s, m in meta.items():
        if m.get("_study") == "TCGA" and m.get("_sample_type") == "Primary Tumor":
            g.setdefault(m.get("primary disease or tissue"), []).append(s)
    return {k: v for k, v in g.items() if k and len(v) >= min_n}

def resolve(symbols):
    """Gene symbol -> versioned Xena ID, via Ensembl REST (the approach used in September)."""
    cache = os.path.join(CACHE, "symbol_to_ensembl.json")
    known = json.load(open(cache)) if os.path.exists(cache) else {}
    todo = [s for s in symbols if s not in known and s not in OVERRIDES]
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    hdr = {"Content-Type": "application/json", "Accept": "application/json"}
    for i in range(0, len(todo), 500):
        batch = todo[i:i + 500]
        for attempt in range(3):
            try:
                r = requests.post(url, headers=hdr, json={"symbols": batch}, timeout=180)
                r.raise_for_status()
                data = r.json()
                break
            except Exception as e:
                print(f"  Ensembl retry {attempt + 1}: {e}")
                time.sleep(10)
        else:
            continue
        for s in batch:
            known[s] = data[s]["id"] if s in data and data[s] else None
    json.dump(known, open(cache, "w"))
    base2ver = load_all_ids()
    out = {}
    for s in symbols:
        b = OVERRIDES.get(s, known.get(s))
        if b and b in base2ver:
            out[s] = base2ver[b]
    return out, [s for s in symbols if s not in out]

def fetch(sym2id, samples, tag, chunk=100):
    """Expression matrix (samples x genes) from Xena, cached to cache/expr_<tag>.pkl."""
    p = os.path.join(CACHE, f"expr_{tag}.pkl")
    if os.path.exists(p):
        df = pd.read_pickle(p)
        if set(sym2id) <= set(df.columns) and set(samples) <= set(df.index):
            return df.loc[samples, list(sym2id)]
    syms, cols = list(sym2id), {}
    for i in range(0, len(syms), chunk):
        part = syms[i:i + chunk]
        for attempt in range(3):
            try:
                res = xena.dataset_probe_values(TOIL, TPM, samples, [sym2id[s] for s in part])
                break
            except Exception as e:
                print(f"  Xena retry {attempt + 1}: {e}")
                time.sleep(15)
        else:
            raise RuntimeError("Xena fetch failed repeatedly")
        for s, vals in zip(part, res[1]):
            cols[s] = pd.to_numeric(pd.Series(vals, index=samples), errors="coerce")
        print(f"  fetched {min(i + chunk, len(syms))}/{len(syms)} genes", flush=True)
    df = pd.DataFrame(cols)
    df.to_pickle(p)
    return df

def bh(p):
    from statsmodels.stats.multitest import multipletests
    p = np.asarray(p, float)
    out = np.full(len(p), np.nan)
    ok = ~np.isnan(p)
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out

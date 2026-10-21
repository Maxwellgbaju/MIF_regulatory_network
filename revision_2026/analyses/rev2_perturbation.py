import re
import pandas as pd
import gseapy as gp

TFS = ["SP1", "RELA", "NFKB1", "HIF1A", "TFAP2A", "EGR1", "KLF4", "IRF1", "IRF2",
       "NFIL3", "ESR1", "TP53", "NR1H2", "RXRA", "MYC"]
LOSS = r"INACTIVAT|KD|KNOCKDOWN|SHRNA|SIRNA|KO|KNOCKOUT|CRISPR|DEPLET|DELET|NULL|LOSS|INHIBIT|SILENC"
GAIN = r"OE|OVEREXPRESS|(?<!IN)ACTIVAT|INDUC|TRANSFECT|GAIN"

libs = gp.get_library_name(organism="Human")
pick = [l for l in libs if re.search(r"Perturbation|CRISPR_KO|shRNA|Knockdown|knockdown", l)]
print("Perturbation libraries found:", pick)

rows, background = [], []
for lib in pick:
    gs = gp.get_library(name=lib, organism="Human")
    up_terms = [t for t in gs if re.search(r"(^|[^A-Z])UP($|[^A-Z])", t.upper())]
    dn_terms = [t for t in gs if re.search(r"(^|[^A-Z])DOWN($|[^A-Z])", t.upper())]
    background.append({"library": lib, "n_terms": len(gs),
                       "MIF_rate_all": sum("MIF" in g for g in gs.values()) / max(len(gs), 1),
                       "MIF_rate_up": sum("MIF" in gs[t] for t in up_terms) / max(len(up_terms), 1),
                       "MIF_rate_down": sum("MIF" in gs[t] for t in dn_terms) / max(len(dn_terms), 1)})
    for term, genes in gs.items():
        T = term.upper()
        for tf in TFS:
            if re.search(rf"(?<![A-Z0-9]){tf}(?![A-Z0-9])", T):
                direction = "UP" if term in up_terms else "DOWN" if term in dn_terms else "unspecified"
                ptype = "loss" if re.search(LOSS, T) else "gain" if re.search(GAIN, T) else "unclear"
                mif = "MIF" in genes
                verdict = ""
                if mif and ptype in ("loss", "gain") and direction in ("UP", "DOWN"):
                    activates = (ptype == "loss" and direction == "DOWN") or (ptype == "gain" and direction == "UP")
                    verdict = "consistent with TF activating MIF" if activates else "consistent with TF repressing MIF"
                species = "mouse/rat" if re.search(r"MOUSE|RAT", T) else "human/unspecified"
                rows.append({"library": lib, "term": term, "TF": tf, "perturbation": ptype,
                             "gene_set_direction": direction, "species": species,
                             "set_size": len(genes), "MIF_in_set": mif, "interpretation": verdict})

hits = pd.DataFrame(rows)
bg = pd.DataFrame(background)
hits.to_csv("results/rev2_perturbation_all_terms.csv", index=False)
bg.to_csv("results/rev2_perturbation_background.csv", index=False)
print("\n=== Background: how often MIF appears in any gene set ===")
print(bg.round(4).to_string(index=False))
if len(hits):
    s = hits.groupby("TF").agg(experiments=("term", "count"), MIF_in_set=("MIF_in_set", "sum")).sort_values("experiments", ascending=False)
    print("\n=== Experiments per TF and how many list MIF ===")
    print(s.to_string())
    print("\n=== Every experiment where MIF responded ===")
    print(hits[hits.MIF_in_set][["TF", "perturbation", "gene_set_direction", "species", "interpretation", "term"]].to_string(index=False))
print("\nSaved results/rev2_*  (absence of MIF does NOT mean no effect: these are top-responding gene sets)")

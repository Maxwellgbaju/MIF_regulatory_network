import pandas as pd

cna = pd.read_csv("cna_6cancers.csv")
muts = pd.read_csv("mutations_6cancers.csv")

study_to_code = {
    "kirc_tcga_pan_can_atlas_2018": "KIRC", "lihc_tcga_pan_can_atlas_2018": "LIHC",
    "stad_tcga_pan_can_atlas_2018": "STAD", "brca_tcga_pan_can_atlas_2018": "BRCA",
    "lusc_tcga_pan_can_atlas_2018": "LUSC", "skcm_tcga_pan_can_atlas_2018": "SKCM",
}
cna["Cancer"] = cna["studyId"].map(study_to_code)
muts["Cancer"] = muts["studyId"].map(study_to_code)

# Total samples per study, derived directly from the CNA fetch (unique sampleIds)
total_samples = cna.groupby("studyId")["sampleId"].nunique().to_dict()
print("Total samples per study (from CNA data):")
for s, n in total_samples.items():
    print(f"  {study_to_code[s]}: {n}")

# --- CNA frequency: % amplified (value>=1) and % deleted (value<=-1) per gene per cancer ---
print("\n=== CNA frequency (%) ===")
cna_summary = []
for study in cna["studyId"].unique():
    n_total = total_samples[study]
    for gene in ["MIF", "SP1", "HIF1A", "RELA"]:
        sub = cna[(cna["studyId"] == study) & (cna["geneSymbol"] == gene)]
        pct_amp = 100 * (sub["value"] >= 1).sum() / n_total
        pct_del = 100 * (sub["value"] <= -1).sum() / n_total
        pct_high_amp = 100 * (sub["value"] == 2).sum() / n_total
        pct_deep_del = 100 * (sub["value"] == -2).sum() / n_total
        cna_summary.append({"Cancer": study_to_code[study], "Gene": gene, "n": n_total,
                             "pct_any_gain": round(pct_amp, 2), "pct_high_amp": round(pct_high_amp, 2),
                             "pct_any_loss": round(pct_del, 2), "pct_deep_del": round(pct_deep_del, 2)})

cna_df = pd.DataFrame(cna_summary)
print(cna_df.to_string(index=False))
cna_df.to_csv("CNA_frequency_summary.csv", index=False)

# --- Mutation frequency: % mutated per gene per cancer ---
print("\n=== Mutation frequency (%) ===")
mut_summary = []
for study in cna["studyId"].unique():  # use same study list to ensure all 6 covered even if 0 mutations
    n_total = total_samples[study]
    for gene in ["MIF", "SP1", "HIF1A", "RELA"]:
        n_mut = muts[(muts["studyId"] == study) & (muts["geneSymbol"] == gene)]["sampleId"].nunique()
        pct_mut = 100 * n_mut / n_total
        mut_summary.append({"Cancer": study_to_code[study], "Gene": gene, "n": n_total,
                             "n_mutated_samples": n_mut, "pct_mutated": round(pct_mut, 2)})

mut_df = pd.DataFrame(mut_summary)
print(mut_df.to_string(index=False))
mut_df.to_csv("Mutation_frequency_summary.csv", index=False)

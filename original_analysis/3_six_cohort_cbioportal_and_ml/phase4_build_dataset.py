import pandas as pd

surv_raw = pd.read_csv("survival_6cancers_raw.csv")
surv = surv_raw.pivot_table(index=["studyId","patientId"], columns="attribute", values="value", aggfunc="first").reset_index()
surv["event"] = surv["OS_STATUS"].astype(str).str.startswith("1").astype(int)
surv["time"] = pd.to_numeric(surv["OS_MONTHS"], errors="coerce")
surv = surv.dropna(subset=["time"])

expr = pd.read_csv("MIF_axis_expression_6cancers.csv")
expr_wide = expr.pivot_table(index=["studyId","sampleId","patientId"], columns="geneSymbol", values="value").reset_index()

merged = expr_wide.merge(surv, on=["studyId","patientId"], how="inner")
kirc = merged[merged["studyId"] == "kirc_tcga_pan_can_atlas_2018"][["patientId","time","event","MIF","SP1","HIF1A","RELA"]].dropna()

print(f"KIRC Phase 4 dataset: n={len(kirc)}, events={kirc['event'].sum()}")
kirc.to_csv("phase4_kirc_dataset.csv", index=False)
print("Saved phase4_kirc_dataset.csv")

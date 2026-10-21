import pandas as pd

cna_freq = pd.read_csv("CNA_frequency_summary.csv")
mut_freq = pd.read_csv("Mutation_frequency_summary.csv")
diploid = pd.read_csv("CNA_diploid_bootstrap.csv")

cancer_order = ["KIRC","LIHC","STAD","BRCA","LUSC","SKCM"]
gene_order = ["MIF","SP1","HIF1A","RELA"]

# --- Table 1: CNA and mutation frequency ---
t1 = cna_freq.merge(mut_freq[["Cancer","Gene","pct_mutated"]], on=["Cancer","Gene"])
t1 = t1[["Cancer","Gene","n","pct_any_gain","pct_high_amp","pct_any_loss","pct_deep_del","pct_mutated"]]
t1.columns = ["Cancer","Gene","n","% any gain","% high amp.","% any loss","% deep del.","% mutated"]
t1["Cancer"] = pd.Categorical(t1["Cancer"], categories=cancer_order, ordered=True)
t1["Gene"] = pd.Categorical(t1["Gene"], categories=gene_order, ordered=True)
t1 = t1.sort_values(["Gene","Cancer"]).reset_index(drop=True)
t1.to_csv("Table_Phase3.4_CNA_mutation_frequency.csv", index=False)

# --- Table 2: diploid-restriction confound test ---
def fmt_p(x):
    return f"{x:.2e}" if x < 0.001 else f"{x:.3f}"

t2 = diploid.copy()
t2.columns = ["Cancer","n (all)","rho (all samples)","n (diploid)","rho (diploid only)",
              "Diff (all-diploid)","Diff 95% CI lower","Diff 95% CI upper","Bootstrap p (shift)"]
for c in ["rho (all samples)","rho (diploid only)","Diff (all-diploid)","Diff 95% CI lower","Diff 95% CI upper"]:
    t2[c] = t2[c].round(3)
t2["Bootstrap p (shift)"] = t2["Bootstrap p (shift)"].apply(fmt_p)
t2["Cancer"] = pd.Categorical(t2["Cancer"], categories=cancer_order, ordered=True)
t2 = t2.sort_values("Cancer").reset_index(drop=True)
t2.to_csv("Table_Phase3.4_diploid_confound_test.csv", index=False)

with pd.ExcelWriter("Phase3.4_Tables.xlsx", engine="openpyxl") as writer:
    t1.to_excel(writer, sheet_name="CNA_Mutation_Frequency", index=False)
    t2.to_excel(writer, sheet_name="Diploid_Confound_Test", index=False)

print("=== Table 1: CNA/mutation frequency ===")
print(t1.to_string(index=False))
print("\n=== Table 2: Diploid-restriction confound test ===")
print(t2.to_string(index=False))
print("\nSaved CSVs and Phase3.4_Tables.xlsx")

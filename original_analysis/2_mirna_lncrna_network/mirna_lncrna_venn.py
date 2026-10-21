import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn2

# All 10 candidate miRNAs target MIF by construction (Step 2.1 selection criterion)
targets_mif = {
    "hsa-miR-654-5p", "hsa-miR-451a", "hsa-miR-1228-5p", "hsa-miR-23a-5p",
    "hsa-miR-363-5p", "hsa-miR-23b-5p", "hsa-miR-541-3p", "hsa-miR-92a-1-5p",
    "hsa-miR-608", "hsa-miR-1537-3p"
}

# miRNAs with ANY lncRNA-sponge evidence (CLIP-seq confirmed OR seed-match)
lncrna_sponge_evidence = {
    "hsa-miR-654-5p", "hsa-miR-541-3p",              # CLIP-seq confirmed (all 4 lncRNAs)
    "hsa-miR-23a-5p", "hsa-miR-363-5p", "hsa-miR-23b-5p",  # seed-match only
    "hsa-miR-1537-3p", "hsa-miR-92a-1-5p"            # seed-match only
}

mif_only = targets_mif - lncrna_sponge_evidence
both = targets_mif & lncrna_sponge_evidence
lncrna_only = lncrna_sponge_evidence - targets_mif  # expected empty

print(f"MIF-only (no lncRNA-sponge evidence): {len(mif_only)} -> {sorted(mif_only)}")
print(f"Both (targets MIF AND has lncRNA-sponge evidence): {len(both)} -> {sorted(both)}")
print(f"lncRNA-only (has lncRNA link but does NOT target MIF): {len(lncrna_only)} -> {sorted(lncrna_only)}")

plt.figure(figsize=(8, 7))
v = venn2([targets_mif, lncrna_sponge_evidence], set_labels=("Targets MIF\n(all 10 candidates)", "lncRNA-sponge evidence\n(CLIP-seq or seed-match)"))
plt.title("Candidate miRNAs: MIF-targeting vs. lncRNA-sponge relationship\n(all miRNAs target MIF by selection criterion; no independent 'lncRNA-only' miRNAs exist in this dataset)", fontsize=10)
plt.tight_layout()
plt.savefig("miRNA_MIF_lncRNA_venn.png", dpi=300, bbox_inches="tight")
plt.savefig("miRNA_MIF_lncRNA_venn.eps", format="eps", bbox_inches="tight")
print("\nSaved miRNA_MIF_lncRNA_venn.png/.eps")

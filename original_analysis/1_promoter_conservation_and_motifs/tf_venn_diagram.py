import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

species_sets = {"human": set(), "chimp": set(), "mouse": set()}
with open("jaspar_scan_results.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        species_sets[row["species"]].add(row["TF_name"])

human, chimp, mouse = species_sets["human"], species_sets["chimp"], species_sets["mouse"]

print(f"Human-only hits: {len(human - chimp - mouse)}")
print(f"Chimp-only hits: {len(chimp - human - mouse)}")
print(f"Mouse-only hits: {len(mouse - human - chimp)}")
print(f"Human & Chimp (not mouse): {len(human & chimp - mouse)}")
print(f"Human & Mouse (not chimp): {len(human & mouse - chimp)}")
print(f"Chimp & Mouse (not human): {len(chimp & mouse - human)}")
print(f"All three (conserved): {len(human & chimp & mouse)}")

plt.figure(figsize=(8, 8))
v = venn3([human, chimp, mouse], set_labels=("Human", "Chimpanzee", "Mouse"))
plt.title("Predicted TF binding sites (JASPAR, 90% threshold)\nby species overlap at the MIF promoter", fontsize=13)
plt.tight_layout()
plt.savefig("TF_venn_3species.png", dpi=300, bbox_inches="tight")
plt.savefig("TF_venn_3species.eps", format="eps", bbox_inches="tight")
print("\nSaved TF_venn_3species.png/.eps")

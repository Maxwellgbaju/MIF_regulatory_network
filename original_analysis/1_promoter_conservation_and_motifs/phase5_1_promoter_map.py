import json
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

TSS = 23894383  # chr22, plus strand
WINDOW_START, WINDOW_END = 23892383, 23896383

# --- Load phyloP conservation track ---
positions, scores = [], []
with open("MIF_conservation_full.bedgraph") as f:
    for line in f:
        chrom, start, end, val = line.strip().split("\t")
        positions.append(int(start))
        scores.append(float(val))
positions = np.array(positions)
scores = np.array(scores)

# --- Load JASPAR hit positions for SP1, HIF1A, RELA (human only, 90% threshold) ---
target_tfs = {"SP1": "#C00000", "HIF1A": "#4472C4", "RELA": "#7030A0"}
tf_positions = {tf: [] for tf in target_tfs}
with open("jaspar_scan_results.csv") as f:
    reader = csv.DictReader(f)
    for row in reader:
        if row["species"] == "human" and row["TF_name"] in target_tfs:
            genomic_pos = WINDOW_START + int(row["position"])
            tf_positions[row["TF_name"]].append(genomic_pos)

# --- Load ENCODE ChIP-seq peaks for SP1 (the only one with confirmed occupancy) ---
data = json.load(open("MIF_ENCODE_TFBS.json"))
sp1_peaks = [(r["chromStart"], r["chromEnd"]) for r in data["encRegTfbsClustered"] if r["name"] == "SP1"]

# --- Build the figure ---
fig, (ax_track, ax_cons) = plt.subplots(2, 1, figsize=(14, 6), gridspec_kw={"height_ratios": [2, 1]}, sharex=True)

# Track panel: TFBS positions as colored ticks, ChIP-seq peaks as thick bars
y_positions = {"SP1": 3, "HIF1A": 2, "RELA": 1}
for tf, color in target_tfs.items():
    y = y_positions[tf]
    for pos in tf_positions[tf]:
        ax_track.plot([pos, pos], [y-0.3, y+0.3], color=color, linewidth=1.2, alpha=0.6)
    ax_track.text(WINDOW_START - 150, y, tf, ha="right", va="center", fontsize=10, fontweight="bold", color=color)

for start, end in sp1_peaks:
    ax_track.add_patch(mpatches.Rectangle((start, y_positions["SP1"]-0.4), end-start, 0.8,
                                            facecolor="#C00000", alpha=0.25, edgecolor="#C00000", linewidth=1.5))

ax_track.axvline(TSS, color="black", linestyle="--", linewidth=1.5)
ax_track.annotate("TSS", xy=(TSS, 3.6), ha="center", fontsize=10, fontweight="bold")
ax_track.set_ylim(0.3, 4)
ax_track.set_yticks([])
ax_track.set_title("MIF Promoter: Transcription Factor Binding Sites (JASPAR predictions, human)\nShaded box = ENCODE ChIP-seq confirmed occupancy (SP1 only)", fontsize=11)

legend_elements = [
    mpatches.Patch(facecolor="#C00000", alpha=0.5, label="SP1 (4-tier evidence: predicted+conserved+TRRUST+ChIP-seq)"),
    mpatches.Patch(facecolor="#4472C4", alpha=0.5, label="HIF1A (3-tier: predicted+conserved+TRRUST)"),
    mpatches.Patch(facecolor="#7030A0", alpha=0.5, label="RELA (3-tier: predicted+conserved+TRRUST)"),
]
ax_track.legend(handles=legend_elements, loc="upper right", fontsize=8)

# Conservation panel: phyloP score across the window
ax_cons.fill_between(positions, scores, 0, where=(scores >= 0), color="#2E7D32", alpha=0.6, step="post")
ax_cons.fill_between(positions, scores, 0, where=(scores < 0), color="#C62828", alpha=0.6, step="post")
ax_cons.axhline(0, color="black", linewidth=0.6)
ax_cons.axvline(TSS, color="black", linestyle="--", linewidth=1.5)
ax_cons.set_ylabel("phyloP\n(conservation)", fontsize=9)
ax_cons.set_xlabel("Genomic position (chr22, hg38)")
ax_cons.set_xlim(WINDOW_START, WINDOW_END)

plt.tight_layout()
plt.savefig("Phase5_1_promoter_map.png", dpi=300, bbox_inches="tight")
plt.savefig("Phase5_1_promoter_map.eps", format="eps", bbox_inches="tight")
print("Saved Phase5_1_promoter_map.png/.eps")

"""Figure 1: study design schematic (no data input)."""
from figstyle import plt, save
from matplotlib.patches import FancyBboxPatch

fig, ax = plt.subplots(figsize=(7.2, 4.8)); ax.axis("off"); ax.set_xlim(0, 100); ax.set_ylim(0, 100)
def box(x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.6,rounding_size=1.5", fc=color, ec="#333333", lw=0.8))
    ax.text(x + w / 2, y + h - 3, title, ha="center", va="top", fontsize=7.4, fontweight="bold")
    ax.text(x + w / 2, y + h - 9, body, ha="center", va="top", fontsize=6.2, linespacing=1.35)
box(1, 52, 30, 45, "1. Candidate nomination",
    "MIF promoter (-2/+2 kb, hg38)\nJASPAR2024 motifs; human,\nchimpanzee and mouse conservation\nTRRUST curation; ENCODE ChIP-seq\nHypothesis-free 12-TF motif screen\nmiRNA screens: multiMiR;\nCyTargetLinker\nlncRNA partners: ENCORI CLIP-seq", "#E8F0FA")
box(35, 52, 30, 45, "2. Expression-based evaluation",
    "Tumour vs normal (TCGA + GTEx)\nCo-expression with MIF:\n32 cancer types, 9,185 tumours\nOverall survival (Cox)\nCross-database machine learning\n(LASSO-Cox, random survival forest)\nInitial six-cohort analysis (cBioPortal)", "#FDF1E2")
box(69, 52, 30, 45, "3. Robustness tests",
    "Three RNA-seq pipelines\n(original TCGA RSEM, PanCanAtlas,\nToil re-quantification)\nTumour-composition adjustment\nTF activity from target genes\nProgramme specificity vs\nexpression-matched genes", "#FBE9E7")
box(18, 3, 30, 39, "4. Direct evidence",
    "ChIP-Atlas: binding within 1 kb\nof the MIF TSS\nPerturbation signatures\n(GEO, CREEDS, LINCS)\nTCGA miRNA-seq test of the\nmiRNA-lncRNA predictions", "#E9F5EC")
box(53, 3, 30, 39, "Integrated evidence",
    "Binding nominates HIF1A and TFAP2A\nas best-supported direct candidates\n\nExpression-based signals explained by\ntumour composition, quantification\nartifacts and generic covariation", "#F2F2F2")
arrow = dict(arrowstyle="-|>", lw=1.1, color="#333")
for x0, x1 in [(31.6, 34.4), (65.6, 68.4)]:
    ax.annotate("", xy=(x1, 75), xytext=(x0, 75), arrowprops=arrow)
ax.annotate("", xy=(33, 42.5), xytext=(16, 51.4), arrowprops=arrow)
ax.annotate("", xy=(68, 42.5), xytext=(84, 51.4), arrowprops=arrow)
ax.annotate("", xy=(52.4, 22), xytext=(48.6, 22), arrowprops=arrow)
save(fig, "Figure1_study_design")

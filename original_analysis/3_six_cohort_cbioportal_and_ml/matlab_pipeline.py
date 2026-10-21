import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Create figure for vector graphics export
fig, ax = plt.subplots(figsize=(8.5, 12), dpi=300)
ax.set_xlim(0, 10)
ax.set_ylim(0, 15)
ax.axis('off')

# Utility to draw rounded box
def draw_box(ax, x, y, width, height, text, bg_color="#FFFFFF", border_color="#4A5568", bold=False):
    box = patches.FancyBboxPatch((x - width/2, y - height/2), width, height,
                                boxstyle="round,pad=0.2,rounding_size=0.15",
                                facecolor=bg_color, edgecolor=border_color, linewidth=1.5)
    ax.add_patch(box)
    weight = 'bold' if bold else 'normal'
    ax.text(x, y, text, ha='center', va='center', fontsize=9, fontfamily='sans-serif', fontweight=weight, color="#2D3748")

def draw_arrow(ax, x1, y1, x2, y2, label=None):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", color="#718096", lw=1.5, mutation_scale=12))
    if label:
        ax.text((x1+x2)/2 + 0.3, (y1+y2)/2, label, fontsize=8, color="#4A5568", fontweight='bold', va='center')

# Outer Cluster Containers
stage1_bg = patches.FancyBboxPatch((1.2, 5.8), 7.6, 7.3, boxstyle="round,pad=0.2,rounding_size=0.1",
                                   facecolor="#F7FAFC", edgecolor="#A0AEC0", linestyle="--", linewidth=1.2)
ax.add_patch(stage1_bg)
ax.text(1.5, 12.8, "Stage I: Conservation-Guided Regulatory Discovery", fontsize=10, fontweight='bold', color="#2D3748")

stage2_bg = patches.FancyBboxPatch((1.2, 1.8), 7.6, 3.4, boxstyle="round,pad=0.2,rounding_size=0.1",
                                   facecolor="#F7FAFC", edgecolor="#A0AEC0", linestyle="--", linewidth=1.2)
ax.add_patch(stage2_bg)
ax.text(1.5, 4.9, "Stage II: Expanded Human-Promoter Screen", fontsize=10, fontweight='bold', color="#2D3748")

# Stage I Nodes
draw_box(ax, 5.0, 14.2, 3.0, 0.5, "MIF Gene (Human)", bg_color="#EDF2F7", bold=True)
draw_box(ax, 5.0, 12.2, 5.5, 0.6, "Human + Chimpanzee + Mouse\nPromoter Conservation Analysis")
draw_box(ax, 5.0, 11.2, 5.5, 0.5, "JASPAR TF Binding Prediction")
draw_box(ax, 5.0, 10.3, 5.5, 0.5, "TRRUST Regulatory Evidence")
draw_box(ax, 5.0, 9.4, 5.5, 0.5, "ENCODE ChIP-seq Validation")
draw_box(ax, 5.0, 8.4, 5.5, 0.6, "High-Confidence Candidates\n(SP1, RELA, HIF1A)", bg_color="#E2E8F0")
draw_box(ax, 5.0, 7.3, 5.5, 0.8, "Clinical Validation (6 TCGA Cohorts)\n• Expression   • Correlation   • Survival\n• Immune Context   • CNA / Mutation")
draw_box(ax, 5.0, 6.2, 5.5, 0.5, "RELA Prioritized (KIRC)", bg_color="#FEB2B2", border_color="#C53030", bold=True)

# Stage II Nodes
draw_box(ax, 5.0, 4.4, 5.5, 0.5, "Human MIF Promoter")
draw_box(ax, 5.0, 3.6, 5.5, 0.5, "JASPAR Ranking (Top TF Candidates)")
draw_box(ax, 5.0, 2.8, 5.5, 0.5, "11 Additional TFs")
draw_box(ax, 5.0, 2.1, 5.5, 0.5, "TCGA-TARGET-GTEx Validation\n• Tumor vs Normal   • Survival   • lncRNA Validation")
draw_box(ax, 5.0, 1.3, 5.5, 0.5, "TFAP2A Prioritized", bg_color="#FEEBC8", border_color="#C05621", bold=True)

# Integrated Modeling Node
draw_box(ax, 5.0, 0.5, 5.5, 0.6, "Integrated Modeling\n• LASSO-Cox   • Stability Selection   • Random Survival Forest")

# Connector Arrows
draw_arrow(ax, 5.0, 13.95, 5.0, 12.5)
draw_arrow(ax, 5.0, 11.9, 5.0, 11.45)
draw_arrow(ax, 5.0, 10.95, 5.0, 10.55)
draw_arrow(ax, 5.0, 10.05, 5.0, 9.65)
draw_arrow(ax, 5.0, 9.15, 5.0, 8.7)
draw_arrow(ax, 5.0, 8.1, 5.0, 7.7)
draw_arrow(ax, 5.0, 6.9, 5.0, 6.45)

draw_arrow(ax, 5.0, 5.95, 5.0, 4.65, label="Next Phase")
draw_arrow(ax, 5.0, 4.15, 5.0, 3.85)
draw_arrow(ax, 5.0, 3.35, 5.0, 3.05)
draw_arrow(ax, 5.0, 2.55, 5.0, 2.35)
draw_arrow(ax, 5.0, 1.85, 5.0, 1.55)
draw_arrow(ax, 5.0, 1.05, 5.0, 0.8)

# Export EPS directly
plt.tight_layout()
plt.savefig("MIF_Pipeline.eps", format='eps', bbox_inches='tight')


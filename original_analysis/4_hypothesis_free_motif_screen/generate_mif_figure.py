import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 1. Define input file name
file_name = "MIF JASPAR DATA_BY_REBCCA.xlsx"

if not os.path.exists(file_name):
    raise FileNotFoundError(f"Could not find '{file_name}' in the current working directory.")

# 2. Load dataset
df = pd.read_excel(file_name)

# 3. Clean and extract the official TF symbol from the JASPAR 'Name' column
df['TF_Symbol'] = df['Name'].apply(lambda x: x.split('.')[-1] if isinstance(x, str) and '.' in x else x)

# 4. Get the Top 12 unique TFs based on maximum JASPAR Score
top_12_df = df.sort_values(by='Score', ascending=False).drop_duplicates(subset=['TF_Symbol']).head(12)
top_12_names = top_12_df['TF_Symbol'].tolist()

# Extract all binding occurrences for these top 12 TFs across the whole promoter sequence
top_12_all_hits = df[df['TF_Symbol'].isin(top_12_names)]

# Determine max promoter sequence length
max_pos = df['End'].max()

# 5. Initialize Figure and Grid Layout (2 Panels)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8.5), gridspec_kw={'height_ratios': [1, 1.2]})
plt.subplots_adjust(hspace=0.35)

# Color Palette (Viridis color map matching both panels)
colors = sns.color_palette("viridis", n_colors=12)

# --- PANEL A: Top 12 TF Ranking Horizontal Bar Chart ---
top_12_sorted = top_12_df.sort_values(by='Score', ascending=True)

bars = ax1.barh(top_12_sorted['TF_Symbol'], top_12_sorted['Score'], color=colors, edgecolor='black', linewidth=0.8)
ax1.set_xlabel('JASPAR Binding Score', fontsize=11, fontweight='bold')
ax1.set_ylabel('Transcription Factor', fontsize=11, fontweight='bold')
ax1.set_title('A. Ranking of Top 12 Candidate Transcription Factors for MIF', fontsize=12, fontweight='bold', loc='left')
ax1.grid(axis='x', linestyle='--', alpha=0.5)

# Label exact scores on the bars
for bar in bars:
    width = bar.get_width()
    ax1.text(width + 0.15, bar.get_y() + bar.get_height()/2, f'{width:.2f}', 
             va='center', ha='left', fontsize=9, fontweight='bold')

ax1.set_xlim(0, max(top_12_sorted['Score']) + 2.5)

# --- PANEL B: Spatial Binding Track Map ---
# Draw promoter base line
ax2.plot([0, max_pos], [0, 0], color='black', linewidth=3, label='MIF Promoter Region')

# Map vertical track index for each TF
tf_y_map = {tf: idx + 1 for idx, tf in enumerate(reversed(top_12_names))}
palette = dict(zip(top_12_names, colors[::-1]))

# Plot binding regions and site midpoints
for tf in top_12_names:
    tf_hits = top_12_all_hits[top_12_all_hits['TF_Symbol'] == tf]
    y_pos = tf_y_map[tf]
    for _, hit in tf_hits.iterrows():
        ax2.plot([hit['Start'], hit['End']], [y_pos, y_pos], color=palette[tf], linewidth=4, solid_capstyle='round')
        ax2.scatter((hit['Start'] + hit['End'])/2, y_pos, color=palette[tf], s=40, zorder=3)

ax2.set_yticks(list(tf_y_map.values()))
ax2.set_yticklabels(list(tf_y_map.keys()), fontsize=10, fontweight='bold')
ax2.set_xlabel('Promoter Position (bp)', fontsize=11, fontweight='bold')
ax2.set_ylabel('Transcription Factor', fontsize=11, fontweight='bold')
ax2.set_title('B. Positional Binding Map Across the MIF Promoter Region', fontsize=12, fontweight='bold', loc='left')
ax2.grid(axis='x', linestyle='--', alpha=0.5)

# Add Transcription Start Site (TSS) dashed line at 2000 bp
if max_pos >= 2000:
    ax2.axvline(2000, color='red', linestyle='--', linewidth=1.5)
    ax2.text(2000, max(tf_y_map.values()) + 0.5, 'TSS', color='red', fontweight='bold', ha='center')

plt.tight_layout()

# Save PNG (Raster, 300 DPI) and EPS (Vector, Journal Standard)
output_png = 'mif_tf_jaspar_figure.png'
output_eps = 'mif_tf_jaspar_figure.eps'

plt.savefig(output_png, dpi=300, bbox_inches='tight')
plt.savefig(output_eps, format='eps', bbox_inches='tight')

print(f"Successfully generated '{output_png}' and '{output_eps}'!")


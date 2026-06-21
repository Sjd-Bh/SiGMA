import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ==============================================================================
# 1. SET PUBLICATION-QUALITY AESTHETICS (Cell Press Style)
# ==============================================================================
# Cell Press generally prefers Arial or Helvetica, vector graphics (PDF/SVG), 
# and clean, uncluttered axes.
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'axes.linewidth': 1.2,
    'xtick.major.width': 1.2,
    'ytick.major.width': 1.2,
    'pdf.fonttype': 42, # Ensures fonts are editable in Illustrator
    'ps.fonttype': 42
})

# ==============================================================================
# 2. GENERATE MOCK DATA (Replace this section with your actual data)
# ==============================================================================
np.random.seed(42)

# --- Mock Data for Breadth of Coverage ---
# Assuming 1 Bulk, 10 MDA, 10 PTA cells
n_cells = 10
breadth_data = pd.DataFrame({
    'Method': ['Bulk'] * 3 + ['MDA'] * n_cells + ['PTA'] * n_cells,
    'Breadth': np.concatenate([
        np.random.normal(95, 1, 3),    # Bulk (high, consistent)
        np.random.normal(55, 10, n_cells), # MDA (lower, highly variable)
        np.random.normal(85, 4, n_cells)   # PTA (higher, more consistent)
    ])
})

# --- Mock Data for Depth of Coverage (Chr 1) ---
positions = np.arange(0, 2500, 1) # e.g., 2500 genomic windows/bins
# Simulate MDA (highly spiky, amplification bias)
mda_depth = np.abs(np.random.normal(10, 15, len(positions)))
# Simulate PTA (flatter, more uniform)
pta_depth = np.abs(np.random.normal(15, 4, len(positions)))

# Simulate centromere gap (coverage drops to 0 in the middle)
centromere_start, centromere_end = 1150, 1350
mda_depth[centromere_start:centromere_end] = 0
pta_depth[centromere_start:centromere_end] = 0

depth_data = pd.DataFrame({
    'Position': positions,
    'MDA_Depth': mda_depth,
    'PTA_Depth': pta_depth
})

# ==============================================================================
# 3. PLOT 1: BREADTH OF COVERAGE (Boxplot + Stripplot)
# ==============================================================================
fig1, ax1 = plt.subplots(figsize=(6, 5))

# Define colors matching the reference images
palette = {'Bulk': '#000000', 'MDA': '#DDA0DD', 'PTA': '#9370DB'} # Adjust colors as needed

# Create the boxplot
sns.boxplot(
    data=breadth_data, x='Method', y='Breadth', 
    ax=ax1, width=0.5, fliersize=0, palette=palette, 
    boxprops=dict(alpha=0.3, edgecolor='black', linewidth=1.5),
    medianprops=dict(color='black', linewidth=1.5),
    whiskerprops=dict(color='black', linewidth=1.5),
    capprops=dict(color='black', linewidth=1.5)
)

# Overlay the individual data points (stripplot)
sns.stripplot(
    data=breadth_data, x='Method', y='Breadth', 
    ax=ax1, size=6, jitter=0.2, palette=palette, edgecolor='black', linewidth=0.5
)

# Calculate Kruskal-Wallis p-value (Comparing MDA and PTA)
group_mda = breadth_data[breadth_data['Method'] == 'MDA']['Breadth']
group_pta = breadth_data[breadth_data['Method'] == 'PTA']['Breadth']
stat, p_val = stats.kruskal(group_mda, group_pta)

# Annotate the p-value on the plot
p_text = f"Kruskal-Wallis (MDA vs PTA), p = {p_val:.1e}" if p_val < 0.001 else f"Kruskal-Wallis, p = {p_val:.3f}"
ax1.text(0.5, 0.95, p_text, transform=ax1.transAxes, ha='center', va='top', fontsize=11)

# Formatting axes
ax1.set_ylabel('Genome breadth (%)', fontweight='bold')
ax1.set_xlabel('')
ax1.set_ylim(0, 105)
sns.despine(ax=ax1, top=True, right=True) # Remove top and right borders

plt.tight_layout()
fig1.savefig('Breadth_of_Coverage.pdf', dpi=300, bbox_inches='tight')
plt.show()

# ==============================================================================
# 4. PLOT 2: MEAN DEPTH OF COVERAGE (Stacked Line Plots)
# ==============================================================================
fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(6, 8), sharex=True, sharey=True)

# Plot MDA (Top Panel) - Using a light blue color like the reference
ax2a.plot(depth_data['Position'], depth_data['MDA_Depth'], color='#5B84C4', linewidth=0.6, alpha=0.9)
ax2a.set_title('MDA', fontweight='bold', fontsize=16)
ax2a.set_ylabel('Mean Depth of Coverage (X)', fontweight='bold')

# Plot PTA (Bottom Panel) - Using a purple color like the reference
ax2b.plot(depth_data['Position'], depth_data['PTA_Depth'], color='#8A5B9C', linewidth=0.6, alpha=0.9)
ax2b.set_title('PTA', fontweight='bold', fontsize=16)
ax2b.set_ylabel('Mean Depth of Coverage (X)', fontweight='bold')
ax2b.set_xlabel('Chr 1 Position', fontweight='bold')

# Formatting for both axes
for ax in [ax2a, ax2b]:
    # Add light gridlines similar to the reference image
    ax.grid(True, axis='both', linestyle='-', color='#E0E0E0', alpha=0.7, zorder=0)
    ax.set_axisbelow(True) # Put grid behind the lines
    
    # Clean up spines (borders) - Reference image uses a light grey bounding box
    for spine in ax.spines.values():
        spine.set_color('#CCCCCC')
        spine.set_linewidth(1.5)
        
    ax.set_ylim(-2, 42) # Adjust based on your actual maximum depth
    
    # Add panel letters (C and D)
    panel_label = 'C' if ax == ax2a else 'D'
    ax.text(-0.15, 1.1, panel_label, transform=ax.transAxes, 
            fontsize=24, fontweight='bold', va='top', ha='right')

plt.tight_layout()
fig2.savefig('Depth_of_Coverage.pdf', dpi=300, bbox_inches='tight')
plt.show()

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats

# ==============================================================================
# 1. SET PUBLICATION-QUALITY AESTHETICS (Cell Press Style)
# ==============================================================================
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
    'pdf.fonttype': 42,
    'ps.fonttype': 42
})

# ==============================================================================
# 2. GENERATE MOCK DATA 
# ==============================================================================
np.random.seed(42)

# --- Mock Data for Breadth of Coverage ---
n_cells = 10
breadth_data = pd.DataFrame({
    'Method': ['Bulk'] * 3 + ['MDA'] * n_cells + ['PTA'] * n_cells,
    'Breadth': np.concatenate([
        np.random.normal(95, 1, 3),    
        np.random.normal(55, 10, n_cells), 
        np.random.normal(85, 4, n_cells)   
    ])
})

# --- Mock Data for Depth of Coverage (Chr 1 - Single Sample) ---
positions = np.arange(0, 2500, 1) 
mda_depth = np.abs(np.random.normal(10, 15, len(positions)))
pta_depth = np.abs(np.random.normal(15, 4, len(positions)))

# Simulate centromere gap 
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

# Define requested colors
#palette = {'Bulk': '#000000', 'MDA': '#DAA520', 'PTA': '#DC143C'} 
colors = {'Bulk': '#000000', 'MDA': '#DAA520', 'PTA': '#DC143C'}

# Create the boxplot with transparent fill and colored edges
sns.boxplot(
    data=breadth_data, x='Method', y='Breadth', 
    ax=ax1, width=0.5, fliersize=0, palette=colors,  
    fill=False, linewidth=1.5  # fill=False makes the box transparent (seaborn >= 0.13)
)

# Overlay the individual data points
sns.stripplot(
    data=breadth_data, x='Method', y='Breadth', legend=False,
    ax=ax1, size=6, jitter=0.2, palette=colors, edgecolor='black', linewidth=0.5
)

# Calculate Kruskal-Wallis p-value 
group_mda = breadth_data[breadth_data['Method'] == 'MDA']['Breadth']
group_pta = breadth_data[breadth_data['Method'] == 'PTA']['Breadth']
stat, p_val = stats.kruskal(group_mda, group_pta)

p_text = f"Kruskal-Wallis (MDA vs PTA), p = {p_val:.1e}" if p_val < 0.001 else f"Kruskal-Wallis, p = {p_val:.3f}"
ax1.text(0.5, 0.95, p_text, transform=ax1.transAxes, ha='center', va='top', fontsize=11)

ax1.set_ylabel('Genome breadth (%)', fontweight='bold')
ax1.set_xlabel('')
ax1.set_ylim(0, 105)
sns.despine(ax=ax1, top=True, right=True) 

plt.tight_layout()
fig1.savefig('Breadth_of_Coverage.pdf', dpi=300, bbox_inches='tight')
plt.show()

# ==============================================================================
# 4. PLOT 2: DEPTH OF COVERAGE FOR A SINGLE SAMPLE (Stacked Line Plots)
# ==============================================================================
fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(6, 8), sharex=True, sharey=True)

# Plot MDA (Top Panel) 
ax2a.plot(depth_data['Position'], depth_data['MDA_Depth'], color='#DAA520', linewidth=0.6, alpha=0.9)
ax2a.set_title('MDA', fontweight='bold', fontsize=16)
ax2a.set_ylabel('Depth of Coverage (X)', fontweight='bold')

# Plot PTA (Bottom Panel) 
ax2b.plot(depth_data['Position'], depth_data['PTA_Depth'], color='#DC143C', linewidth=0.6, alpha=0.9)
ax2b.set_title('PTA', fontweight='bold', fontsize=16)
ax2b.set_ylabel('Depth of Coverage (X)', fontweight='bold')
ax2b.set_xlabel('Chr 1 Position', fontweight='bold')

# Formatting for both axes
for ax in [ax2a, ax2b]:
    ax.grid(True, axis='both', linestyle='-', color='#E0E0E0', alpha=0.7, zorder=0)
    ax.set_axisbelow(True) 
    
    for spine in ax.spines.values():
        spine.set_color('#CCCCCC')
        spine.set_linewidth(1.5)
        
    ax.set_ylim(-2, 42) 
    
    panel_label = 'C' if ax == ax2a else 'D'
    ax.text(-0.15, 1.1, panel_label, transform=ax.transAxes, 
            fontsize=24, fontweight='bold', va='top', ha='right')

plt.tight_layout()
fig2.savefig('Depth_of_Coverage.pdf', dpi=300, bbox_inches='tight')
plt.show()

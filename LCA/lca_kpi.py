"""
lca_kpi.py
==========
Grouped bar chart of LCA KPI results for IL and Water sorbent systems.
Multiplicity auto-calculated from raw values so all KPIs plot on same scale.
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

KPI_XLSX   = r'LCAResults_KPIs.xlsx'
OUTPUT_DIR = r'.'
DPI        = 200

WATER_COLOR = '#2166AC'
IL_COLOR    = '#E48B17'

plt.rcParams.update({
    'font.family':      'sans-serif',
    'font.size':        13,
    'axes.titlesize':   14,
    'axes.labelsize':   13,
    'xtick.labelsize':  11,
    'ytick.labelsize':  12,
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
})

# Load data
df = pd.read_excel(KPI_XLSX, header=0)
df.columns = df.columns.str.strip()

n_kpi = len(df)
bar_w = 0.32
x     = np.arange(n_kpi)

fig, ax = plt.subplots(figsize=(16, 7))
fig.patch.set_facecolor('white')

il_scaled    = []
water_scaled = []
x_labels     = []

for _, row in df.iterrows():
    il_raw    = float(row['IL'])
    water_raw = float(row['Water'])

    # Auto-calculate multiplicity from max of IL and water
    max_val = max(abs(il_raw), abs(water_raw))
    if max_val == 0:
        exp = 0
    else:
        exp = -int(np.floor(np.log10(max_val)))

    mult = 10 ** exp
    il_scaled.append(il_raw * mult)
    water_scaled.append(water_raw * mult)

    # Format multiplicity for label
    if exp == 0:
        mult_str = ''
    elif exp > 0:
        mult_str = f' ($\\times 10^{{-{exp}}}$)'
    else:
        mult_str = f' ($\\times 10^{{{abs(exp)}}}$)'

    unit_str = str(row['Units']).strip()
    unit_str = str(row['Units']).strip()
    if mult_str:
        label = f"{row['Impact Category']}\n{unit_str}\n{mult_str}"
    else:
        label = f"{row['Impact Category']}\n{unit_str}"
    x_labels.append(label)

il_scaled    = np.array(il_scaled)
water_scaled = np.array(water_scaled)

# Plot paired bars — Water left, IL right
ax.bar(x - bar_w/2, water_scaled, bar_w,
       color=WATER_COLOR, hatch='',
       edgecolor='white', linewidth=0.5,
       alpha=0.85, zorder=3, label='Water')
ax.bar(x + bar_w/2, il_scaled, bar_w,
       color=IL_COLOR, hatch='//',
       edgecolor='white', linewidth=0.5,
       alpha=0.85, zorder=3, label='IL ([bmim][DCA])')

# Formatting
ax.set_xticks(x)
ax.set_xticklabels(x_labels, ha='right', linespacing=1.3,
                   rotation=30, multialignment='center')
ax.tick_params(axis='x', pad=8)
ax.set_ylabel('Scaled Impact Value', fontsize=13)
ax.yaxis.grid(True, color='#CCCCCC', linestyle='--',
              linewidth=0.5, alpha=0.7, zorder=0)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlim(-0.6, n_kpi - 0.4)
ax.set_ylim(bottom=0)

# Legend
wat_patch = mpatches.Patch(facecolor=WATER_COLOR, hatch='',
                            edgecolor='#444444', linewidth=0.8,
                            label='Water Sorbent')
il_patch  = mpatches.Patch(facecolor=IL_COLOR, hatch='//',
                            edgecolor='#444444', linewidth=0.8,
                            label='IL Sorbent ([bmim][DCA])')
ax.legend(handles=[wat_patch, il_patch], fontsize=12,
          loc='upper right', frameon=True,
          framealpha=0.95, edgecolor='#cccccc')

plt.tight_layout()
plt.subplots_adjust(bottom=0.35)
fname = os.path.join(OUTPUT_DIR, 'LCA_KPI.png')
plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {fname}')

"""
lca_tornado.py
==============
Paired bar tornado plot for LCA sensitivity analysis.
IL = orange, Water = blue. Linear scale on x-axis.
NaN values (sensitivity not applicable) are skipped per sorbent.

Usage: python lca_tornado.py
Requirements: pip install pandas openpyxl matplotlib numpy
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

TORNADO_XLSX = r'LCAResults_TornadoPlot.xlsx'
OUTPUT_DIR   = r'.'
DPI          = 200

WATER_COLOR = '#2166AC'
IL_COLOR    = '#E48B17'
BASELINE_COLOR = '#333333'

plt.rcParams.update({
    'font.family':      'sans-serif',
    'font.size':        11,
    'axes.titlesize':   12,
    'axes.labelsize':   11,
    'xtick.labelsize':  9,
    'ytick.labelsize':  10,
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
})

# Load data
df = pd.read_excel(TORNADO_XLSX, header=0)
df.columns = df.columns.str.strip()

# Get IL and Water subsets
il_df  = df[df['Sorbent'] == 'IL'].copy().reset_index(drop=True)
wat_df = df[df['Sorbent'] == 'Water'].copy().reset_index(drop=True)

# Get baselines
il_base  = il_df['Baseline GWP'].iloc[0]
wat_base = wat_df['Baseline GWP'].iloc[0]

# Compute deltas
for d in [il_df, wat_df]:
    base = d['Baseline GWP'].iloc[0]
    d['delta_low']  = (d['Low GWP'].astype(float)  - base) / abs(base) * 100
    d['delta_high'] = (d['High GWP'].astype(float) - base) / abs(base) * 100
    d['max_range']  = d[['delta_low','delta_high']].abs().max(axis=1)

# Sort by IL max range
il_df = il_df.sort_values('max_range', ascending=True).reset_index(drop=True)
var_order = il_df['Variable'].tolist()

# Align water to same order
wat_rows = []
for v in var_order:
    match = wat_df[wat_df['Variable'] == v]
    wat_rows.append(match.iloc[0] if len(match) > 0 else None)

n_vars = len(var_order)
y_pos  = np.arange(n_vars)
bar_h  = 0.35

fig, ax = plt.subplots(figsize=(11, max(5, n_vars * 1.1)))

for i in range(n_vars):
    il_row = il_df.iloc[i]
    pairs  = [(il_row, IL_COLOR, bar_h/2, '//')]
    if wat_rows[i] is not None:
        pairs.append((wat_rows[i], WATER_COLOR, -bar_h/2, ''))

    for row, color, yoff, hatch in pairs:
        d_low  = row['delta_low']
        d_high = row['delta_high']

        if pd.notna(d_high) and d_high > 0:
            ax.barh(y_pos[i] + yoff, d_high, bar_h,
                    left=0, color=color, hatch='//',
                    alpha=0.85, edgecolor='white',
                    linewidth=0.5, zorder=3)
            ax.text(d_high + 0.5, y_pos[i] + yoff,
                    f'{d_high:+.1f}%',
                    va='center', ha='left', fontsize=8, color='black')

        if pd.notna(d_low) and d_low < 0:
            ax.barh(y_pos[i] + yoff, d_low, bar_h,
                    left=0, color=color, hatch='',
                    alpha=0.85, edgecolor='white',
                    linewidth=0.5, zorder=3)
            ax.text(d_low - 0.5, y_pos[i] + yoff,
                    f'{d_low:+.1f}%',
                    va='center', ha='right', fontsize=8, color='black')

# Baseline lines
ax.axvline(0, color=BASELINE_COLOR, linewidth=1.5,
           linestyle='--', alpha=0.8, zorder=4)

from matplotlib.offsetbox import TextArea, HPacker, VPacker, AnnotationBbox
from matplotlib.transforms import blended_transform_factory

def _fmt(x):
    if pd.isna(x):
        return 'N/A'
    return str(x)

ax.set_yticks(y_pos)
ax.set_yticklabels([])
ax.tick_params(axis='y', which='major', length=4)

trans = blended_transform_factory(ax.transAxes, ax.transData)

for i, (_, row) in enumerate(il_df.iterrows()):
    y   = y_pos[i]
    var = row['Variable']
    lo  = _fmt(row['Variable Low'])
    bl  = _fmt(row['Variable Baseline'])
    hi  = _fmt(row['Variable High'])

    t_var = TextArea(var, textprops=dict(size=10))
    t_lo  = TextArea(f"({lo}, ", textprops=dict(size=8, color='#444444'))
    if hi == 'N/A':
        t_bl  = TextArea(f"{bl})", textprops=dict(size=8, color='#444444', fontweight='bold'))
        line2 = HPacker(children=[t_lo, t_bl], pad=0, sep=0)
    else:
        t_bl  = TextArea(bl,          textprops=dict(size=8, color='#444444', fontweight='bold'))
        t_hi  = TextArea(f", {hi})",  textprops=dict(size=8, color='#444444'))
        line2 = HPacker(children=[t_lo, t_bl, t_hi], pad=0, sep=0)
    vbox  = VPacker(children=[t_var, line2],      pad=0, sep=2, align='right')

    ab = AnnotationBbox(vbox, xy=(0, y), xycoords=trans,
                        xybox=(-8, 0), boxcoords='offset points',
                        box_alignment=(1.0, 0.5),
                        frameon=False, zorder=10)
    ab.set_clip_on(False)
    ax.add_artist(ab)

xmax = max(
    (abs(df['High GWP'].fillna(df['Baseline GWP']) - df['Baseline GWP']) / df['Baseline GWP'].abs()).max(),
    (abs(df['Low GWP'] - df['Baseline GWP']) / df['Baseline GWP'].abs()).max()
) * 100 * 1.3
# leave room left of the largest decrease for its % label
xmin = (abs(df['Low GWP'] - df['Baseline GWP']) / df['Baseline GWP'].abs()).max() * 100 * 1.3
ax.set_xlim(-max(100, xmin), xmax)
ax.set_xlabel('% Change in GWP from Baseline', fontsize=10)

ax.xaxis.grid(True, color='#CCCCCC', linestyle='--',
              linewidth=0.5, alpha=0.7, zorder=0)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
il_patch  = mpatches.Patch(facecolor=IL_COLOR, hatch='',
                            edgecolor='#444444', linewidth=0.8,
                            label=f'IL ([bmim][DCA]) — baseline {il_base:.1f} kg CO₂-eq/kg CO₂')
wat_patch = mpatches.Patch(facecolor=WATER_COLOR, hatch='',
                            edgecolor='#444444', linewidth=0.8,
                            label=f'Water — baseline {wat_base:.1f} kg CO₂-eq/kg CO₂')
inc_patch = mpatches.Patch(facecolor='grey', hatch='//',
                            edgecolor='#444444', linewidth=0.8,
                            alpha=0.7, label='Increases GWP (hatched)')
dec_patch = mpatches.Patch(facecolor='grey', hatch='',
                            edgecolor='#444444', linewidth=0.8,
                            alpha=0.7, label='Decreases GWP (solid)')
ax.legend(handles=[il_patch, wat_patch, inc_patch, dec_patch],
          fontsize=9, loc='lower right',
          frameon=True, framealpha=0.95, edgecolor='#cccccc')

plt.tight_layout()
fname = os.path.join(OUTPUT_DIR, 'LCA_tornado.png')
plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {fname}')

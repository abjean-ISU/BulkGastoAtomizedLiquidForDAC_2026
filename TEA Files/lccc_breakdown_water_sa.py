"""
lccc_breakdown_water_sa.py
==========================
LCCC breakdown for 3 water sorbent sensitivity cases:
  Water Baseline | Low Liquid Holdup | High Liquid Holdup

Left panel : absolute LCCC ($/kg CO2), continuous y-axis
Right panel: 100% contribution breakdown

Colors are from the Okabe-Ito palette (colorblind-safe), blue/green family.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

plt.rcParams.update({
    'font.family':     'sans-serif',
    'font.size':       13,
    'axes.labelsize':  13,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
})

# ============================================================
# USER SETTINGS
# ============================================================
XLSX_PATH  = r'C:\Users\anbarron\Box\Research\Projects\DAC\Modeling\TEA\LCCCBreakdown_WaterLHSA.xlsx'
OUTPUT_DIR = r'C:\Users\anbarron\Box\Research\Projects\DAC\Modeling\TEA'
DPI        = 200

# ============================================================
# COLORS — Okabe-Ito blue/teal/green, all colorblind-safe
# ============================================================
BASE_COLORS = {
    'Water Baseline': '#0072B2',   # Okabe-Ito blue
    'Low LH':         '#009E73',   # Okabe-Ito bluish green
    'High LH':        '#4DAF4A',   # ColorBrewer Set1 green — brighter, clearly distinct
}
CASE_LABELS = {
    'Water Baseline': 'Baseline (0.4% Liquid Holdup)',
    'Low LH':         '0.04% Liquid Holdup',
    'High LH':        '8% Liquid Holdup',
}
HATCHES = ['', '//', '\\\\', 'xx', '..', '++', '--']

def make_shades(base_hex, n):
    import matplotlib.colors as mc
    rgb = mc.to_rgb(base_hex)
    return [tuple(1 - (1 - c) * f for c in rgb)
            for f in np.linspace(0.40, 0.92, n)]

# ============================================================
# LOAD & ORDER DATA
# ============================================================
df = pd.read_excel(XLSX_PATH, header=0, index_col=0)
df = df[df.index.notna()]   # drop unlabeled total/fraction rows below the table
df.columns = ['Water Baseline', 'Low LH', 'High LH']

desired_order = [
    'Annualized Capital',
    'LO+M',
    'Cooling Water',
    'Electricity',
    'Property Tax + Insurance',
    'Refrigeration',
    'Raw Materials (Sorbent)',
]
df         = df.reindex(desired_order)
components = df.index.tolist()
n_comp     = len(components)
cases      = list(df.columns)

vals   = {c: df[c].values.astype(float) for c in cases}
shades = {c: make_shades(BASE_COLORS[c], n_comp) for c in cases}
totals = {c: float(df[c].sum()) for c in cases}

# ============================================================
# FIGURE — (a) absolute LCCC   (b) 100% contribution
# ============================================================
bar_width = 0.50
x = np.array([0, 1, 2])

fig, (ax_left, ax_pct) = plt.subplots(1, 2, figsize=(15, 8))
fig.patch.set_facecolor('white')
fig.subplots_adjust(wspace=0.35, bottom=0.28)

for ax in (ax_left, ax_pct):
    ax.set_facecolor('white')

# ============================================================
# LEFT PANEL — absolute LCCC, single continuous axis
# ============================================================
bottoms = {c: 0.0 for c in cases}
legend_handles = []

for i, comp in enumerate(components):
    hatch = HATCHES[i % len(HATCHES)]
    for j, case in enumerate(cases):
        ax_left.bar(x[j], vals[case][i], bar_width,
                    bottom=bottoms[case],
                    color=shades[case][i],
                    hatch=hatch,
                    edgecolor='white',
                    linewidth=0.5,
                    zorder=3)
        bottoms[case] += vals[case][i]

    legend_handles.append(mpatches.Patch(
        facecolor='#888888',
        hatch=hatch,
        edgecolor='#444444',
        linewidth=0.8,
        label=comp,
    ))

max_total = max(totals.values())
ax_left.set_ylim(0, max_total * 1.12)
ax_left.set_xlim(-0.45, 2.45)

# Total annotations
for j, case in enumerate(cases):
    tot = totals[case]
    ax_left.text(x[j], tot + max_total * 0.01,
                 f'${tot:.0f}/kg CO\u2082',
                 ha='center', va='bottom', fontsize=11,
                 fontweight='bold', color='black')

ax_left.set_xticks(x)
ax_left.set_xticklabels(
    ['Baseline\n(0.4% Liquid Holdup)', '0.04% Liquid\nHoldup', '8% Liquid\nHoldup'],
    fontsize=12)
ax_left.set_ylabel('LCCC ($/kg CO\u2082 captured)', fontsize=13)
ax_left.yaxis.grid(True, color='#CCCCCC', linestyle='--',
                   linewidth=0.5, alpha=0.7, zorder=0)
ax_left.set_axisbelow(True)
ax_left.spines['top'].set_visible(False)
ax_left.spines['right'].set_visible(False)

# ============================================================
# RIGHT PANEL — 100% contribution
# ============================================================
pct_bottoms = {c: 0.0 for c in cases}
for i, comp in enumerate(components):
    hatch = HATCHES[i % len(HATCHES)]
    for j, case in enumerate(cases):
        pct = vals[case][i] / totals[case] * 100
        ax_pct.bar(x[j], pct, bar_width,
                   bottom=pct_bottoms[case],
                   color=shades[case][i],
                   hatch=hatch,
                   edgecolor='white',
                   linewidth=0.5,
                   zorder=3)
        pct_bottoms[case] += pct

ax_pct.set_ylim(0, 100)
ax_pct.set_xlim(-0.45, 2.45)
ax_pct.set_xticks(x)
ax_pct.set_xticklabels(
    ['Baseline\n(0.4% Liquid Holdup)', '0.04% Liquid\nHoldup', '8% Liquid\nHoldup'],
    fontsize=12)
ax_pct.set_ylabel('Contribution to LCCC (%)', fontsize=13)
ax_pct.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
ax_pct.yaxis.grid(True, color='#CCCCCC', linestyle='--',
                  linewidth=0.5, alpha=0.7, zorder=0)
ax_pct.set_axisbelow(True)
ax_pct.spines['top'].set_visible(False)
ax_pct.spines['right'].set_visible(False)

# ============================================================
# LEGEND — placed below both panels with clear gap
# ============================================================
color_keys = [
    mpatches.Patch(facecolor=BASE_COLORS[c], edgecolor='white',
                   label=CASE_LABELS[c])
    for c in cases
]
fig.legend(
    handles=legend_handles + color_keys,
    title='Cost Component',
    title_fontsize=12,
    fontsize=11,
    loc='lower center',
    bbox_to_anchor=(0.5, 0.0),
    ncol=3,
    frameon=True,
    framealpha=0.95,
    edgecolor='#cccccc',
    handlelength=3.5,
    handleheight=2.0,
)

# ============================================================
# PANEL LABELS
# ============================================================
pos_left = ax_left.get_position()
pos_pct  = ax_pct.get_position()
label_y  = max(pos_left.y1, pos_pct.y1) + 0.02
fig.text(pos_left.x0 - 0.01, label_y, '(a)',
         fontsize=14, fontweight='bold', va='bottom', ha='right',
         transform=fig.transFigure)
fig.text(pos_pct.x0 - 0.01, label_y, '(b)',
         fontsize=14, fontweight='bold', va='bottom', ha='right',
         transform=fig.transFigure)

# ============================================================
# SAVE
# ============================================================
fname = os.path.join(OUTPUT_DIR, 'LCCC_breakdown_water_sa.png')
plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {fname}')

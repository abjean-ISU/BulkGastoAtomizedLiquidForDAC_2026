"""
combined_tornado.py — Paired bar tornado with symlog percent scale
IL = orange, Water = blue, textures = cost components in legend
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

IL_XLSX    = r'ILTornadoPlotValues.xlsx'
WAT_XLSX   = r'WaterTornadoPlotValues.xlsx'
OUTPUT_DIR = r'.'
DPI        = 200

WATER_COLOR = '#2166AC'
IL_COLOR    = '#E48B17'
BASELINE_COLOR = '#333333'

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 13,
    'axes.titlesize': 14, 'axes.labelsize': 13,
    'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
})

def load_tornado(path, baseline_override=None):
    df = pd.read_excel(path, header=0)
    df.columns = ['Variable','Low','Baseline','High',
                  'Low_LCCC','Baseline_LCCC','High_LCCC']
    df['Variable'] = df['Variable'].str.strip()
    df = df.dropna(subset=['Low_LCCC','High_LCCC'], how='all').reset_index(drop=True)
    baseline = baseline_override if baseline_override else float(df['Baseline_LCCC'].iloc[0])
    df['delta_low']  = (df['Low_LCCC'].astype(float)  - baseline) / baseline * 100
    df['delta_high'] = (df['High_LCCC'].astype(float) - baseline) / baseline * 100
    df['delta_low']  = df['delta_low'].fillna(0)
    df['delta_high'] = df['delta_high'].fillna(0)
    return df, baseline

il_df,  il_base  = load_tornado(IL_XLSX)
wat_df, wat_base = load_tornado(WAT_XLSX)

# Remove liquid holdup from both datasets
EXCLUDE_VARS = ['Liquid Holdup']
il_df  = il_df[~il_df['Variable'].isin(EXCLUDE_VARS)].reset_index(drop=True)
wat_df = wat_df[~wat_df['Variable'].isin(EXCLUDE_VARS)].reset_index(drop=True)

def fmt_val(v):
    """Format a parameter value: numbers to 3 sig figs, strings as-is."""
    s = str(v).strip()
    try:
        return f'{float(s):.3g}'
    except ValueError:
        return s

def fmt_baseline_bold(s):
    """Return baseline value formatted bold for matplotlib tick labels.
    Moves % outside of math mode to avoid mathtext parse errors."""
    if '%' in s:
        idx = s.find('%')
        before, after = s[:idx], s[idx:]
        return f'$\\mathbf{{{before}}}${after}' if before else s
    return f'$\\mathbf{{{s}}}$'

def escape_mpl(s):
    """Escape $ so it renders as a literal dollar sign outside math mode."""
    return s.replace('$', r'\$')

# Shared variables only (water variables are subset of IL)
shared_vars = wat_df['Variable'].tolist()

# Sort all IL variables by range
il_df['max_range'] = il_df[['delta_low','delta_high']].abs().max(axis=1)
il_df = il_df.sort_values('max_range', ascending=True).reset_index(drop=True)
var_order = il_df['Variable'].tolist()

# Water only has data for shared variables — None for IL-only
wat_rows = []
for v in var_order:
    match = wat_df[wat_df['Variable'] == v]
    wat_rows.append(match.iloc[0] if len(match) > 0 else None)

n_vars = len(var_order)
y_pos  = np.arange(n_vars)
bar_h  = 0.35

fig, ax = plt.subplots(figsize=(11, max(5, n_vars * 0.9)))

HALF = bar_h / 2

for i in range(n_vars):
    pairs = [(il_df.iloc[i], IL_COLOR, HALF/2)]
    if wat_rows[i] is not None:
        pairs.append((wat_rows[i], WATER_COLOR, -HALF/2))
    for row, color, yoff in pairs:
        d_low  = row['delta_low']
        d_high = row['delta_high']
        y = y_pos[i] + yoff

        # Low perturbation — solid
        if d_low != 0:
            ax.barh(y, d_low, HALF,
                    left=0, color=color, alpha=0.85,
                    edgecolor='white', linewidth=0.5, zorder=3)
            ha  = 'left' if d_low > 0 else 'right'
            off = max(abs(d_low) * 0.05, 1.0) * (1 if d_low > 0 else -1)
            xlab = d_low + off
            # both ends on the same side: put the shorter bar's label on the other side of zero
            if d_low * d_high > 0 and abs(d_low) < abs(d_high):
                ha   = 'right' if d_low > 0 else 'left'
                xlab = -1.0 if d_low > 0 else 1.0
            ax.text(xlab, y, f'{d_low:+.1f}%',
                    va='center', ha=ha, fontsize=11, color='black')

        # High perturbation — hatched
        if d_high != 0:
            ax.barh(y, d_high, HALF,
                    left=0, color=color, alpha=0.85,
                    hatch='//', edgecolor='white', linewidth=0.5, zorder=4)
            ha  = 'left' if d_high > 0 else 'right'
            off = max(abs(d_high) * 0.05, 1.0) * (1 if d_high > 0 else -1)
            xlab = d_high + off
            if d_low * d_high > 0 and abs(d_high) < abs(d_low):
                ha   = 'right' if d_high > 0 else 'left'
                xlab = -1.0 if d_high > 0 else 1.0
            ax.text(xlab, y, f'{d_high:+.1f}%',
                    va='center', ha=ha, fontsize=11, color='black')

ax.axvline(0, color=BASELINE_COLOR, linewidth=1.5,
           linestyle='--', alpha=0.8, zorder=4)

yticklabels = []
for v in var_order:
    row = il_df[il_df['Variable'] == v].iloc[0]
    low_s  = fmt_val(row['Low'])
    base_s = fmt_val(row['Baseline'])
    high_s = fmt_val(row['High'])
    v_lbl    = escape_mpl(v)
    low_lbl  = escape_mpl(low_s)
    high_lbl = escape_mpl(high_s)
    yticklabels.append(f"{v_lbl}\n({low_lbl}, {fmt_baseline_bold(base_s)}, {high_lbl})")

ax.set_yticks(y_pos)
ax.set_yticklabels(yticklabels)

ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:+.0f}%' if x != 0 else '0%'))

ax.set_xlabel('% Change in LCCC ($/kg CO₂) from baseline', fontsize=13)
ax.xaxis.grid(True, color='#CCCCCC', linestyle='--',
              linewidth=0.5, alpha=0.7, zorder=0)
ax.set_axisbelow(True)
ax.set_xscale('symlog', linthresh=10)
ax.set_xlim(left=-1000)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Legend
il_patch  = mpatches.Patch(color=IL_COLOR,    alpha=0.85, label=f'IL ([bmim][DCA]) — baseline ${il_base:.0f}/kg CO₂')
wat_patch = mpatches.Patch(color=WATER_COLOR, alpha=0.85, label=f'Water — baseline ${wat_base:.0f}/kg CO₂')
low_patch  = mpatches.Patch(facecolor='grey', hatch='',   edgecolor='#444444',
                            alpha=0.7, label='Solid — Low perturbation')
high_patch = mpatches.Patch(facecolor='grey', hatch='//', edgecolor='#444444',
                             alpha=0.7, label='Hatched — High perturbation')
ax.legend(handles=[il_patch, wat_patch, low_patch, high_patch],
          fontsize=11, loc='upper center', bbox_to_anchor=(0.5, -0.14), ncol=2, framealpha=0.9, edgecolor='#cccccc')

plt.tight_layout()
fname = os.path.join(OUTPUT_DIR, 'Combined_tornado.png')
plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {fname}')

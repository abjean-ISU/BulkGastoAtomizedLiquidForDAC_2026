"""
il_surface.py
=============
Three-panel heatmap with contour overlay showing LCCC as a function
of IL unit cost and sorbent loss rate across three absorption efficiency
scenarios (Baseline, High, Low).

Usage: python il_surface.py
Requirements: pip install pandas openpyxl matplotlib numpy scipy
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import os

# ============================================================
# USER SETTINGS
# ============================================================
XLSX_PATH  = r'ILSobentLossCostSurfacePlot.xlsx'
OUTPUT_DIR = r'.'
DPI        = 200

BASELINE_IL_COST   = 30      # $/kg — baseline marker
BASELINE_LOSS_RATE = 0.005   # fraction

CONTOUR_LEVELS = [100, 500, 1000, 2000, 4000, 7000, 13000]

SHEET_TITLES = {
    'Low Abs Eff':      'Low Absorption Capacity',
    'Baseline Abs Eff': 'Baseline Absorption Capacity',
    'High Abs Eff':     'High Absorption Capacity',
}

plt.rcParams.update({
    'font.family':      'sans-serif',
    'font.size':        13,
    'axes.titlesize':   13,
    'axes.labelsize':   13,
    'xtick.labelsize':  12,
    'ytick.labelsize':  12,
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
})

# ============================================================
# LOAD DATA FROM ALL THREE SHEETS
# ============================================================
def load_sheet(path, sheet):
    df = pd.read_excel(path, sheet_name=sheet, header=0, index_col=0)
    # Parse IL cost column headers
    il_costs = []
    for col in df.columns:
        c = str(col).replace('$','').replace('/kg','').strip()
        il_costs.append(float(c))
    # Parse loss rate row indices
    loss_rates = []
    for idx in df.index:
        loss_rates.append(float(str(idx).strip()))
    Z = df.values.astype(float)
    return np.array(il_costs), np.array(loss_rates), Z

xl = pd.ExcelFile(XLSX_PATH)
sheets = ['Low Abs Eff', 'Baseline Abs Eff', 'High Abs Eff']
# Compute global min/max for consistent colorscale across panels
all_Z = []
for sheet in sheets:
    _, _, Z = load_sheet(XLSX_PATH, sheet)
    all_Z.append(Z)
z_min = min(Z.min() for Z in all_Z)
z_max = max(Z.max() for Z in all_Z)

# ============================================================
# FIGURE — 1 row x 3 panels
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5.5),
                          sharey=True, sharex=True)

norm = mcolors.LogNorm(vmin=max(z_min, 1), vmax=z_max)

for ax, sheet in zip(axes, sheets):
    X_vals, Y_vals, Z = load_sheet(XLSX_PATH, sheet)

    # Interpolate onto fine grid
    interp = RegularGridInterpolator(
        (Y_vals, X_vals), Z,
        method='linear', bounds_error=False, fill_value=None
    )
    X_fine = np.linspace(X_vals.min(), X_vals.max(), 200)
    Y_fine = np.logspace(np.log10(Y_vals.min()),
                          np.log10(Y_vals.max()), 200)
    Xg, Yg = np.meshgrid(X_fine, Y_fine)
    Zg = interp((Yg, Xg))

    # Filled heatmap
    pcm = ax.pcolormesh(Xg, Yg * 100, Zg,
                        cmap='YlOrBr', norm=norm,
                        shading='auto', zorder=1)

    # Contour lines
    levels = [l for l in CONTOUR_LEVELS
              if z_min <= l <= z_max]
    if levels:
        cs = ax.contour(Xg, Yg * 100, Zg,
                        levels=levels,
                        colors='black',
                        linewidths=1.0,
                        alpha=0.85, zorder=2)
        ax.clabel(cs, fmt='$%.0f', fontsize=10,
                  inline=True, inline_spacing=4)

    # Baseline marker
    ax.plot(BASELINE_IL_COST, BASELINE_LOSS_RATE * 100,
            marker='*', markersize=14,
            color='#2166AC', markeredgecolor='white',
            markeredgewidth=1.0, zorder=5,
            label=f'Baseline (${BASELINE_IL_COST}/kg, '
                  f'{BASELINE_LOSS_RATE*100:.1f}%)')

    ax.set_title(SHEET_TITLES[sheet], fontsize=13)

    # Log scale y-axis
    ax.set_yscale('log')
    ax.set_yticks([y * 100 for y in Y_vals])
    ax.set_yticklabels([f'{y*100:.3f}%' if y*100 < 0.01
                        else f'{y*100:.2f}%' for y in Y_vals])
    ax.yaxis.set_minor_formatter(plt.NullFormatter())

    # X ticks
    ax.set_xticks(X_vals)
    ax.set_xticklabels([f'${int(v)}' for v in X_vals])

    ax.set_xlabel('IL Unit Cost ($/kg)', fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    if ax == axes[0]:
        ax.set_ylabel('Sorbent Loss Rate (%)', fontsize=13)
        ax.legend(fontsize=11, loc='upper left',
                  framealpha=0.9, edgecolor='#cccccc')

# Shared colorbar on the right
fig.subplots_adjust(right=0.88)
cbar = fig.colorbar(pcm, ax=axes.tolist(), pad=0.02, fraction=0.015,
                    location='right')
cbar.set_label('LCCC ($/kg CO₂ captured)', fontsize=13)
cbar.ax.tick_params(labelsize=12)
cbar.ax.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: f'${x:,.0f}'))

fname = os.path.join(OUTPUT_DIR, 'IL_surface.png')
plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {fname}')

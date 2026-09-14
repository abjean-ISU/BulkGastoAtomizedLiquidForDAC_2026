"""
combined_foak_noak.py
=====================
FOAK vs NOAK — Water always left of IL within each group.
IL = orange shades, Water = blue shades.
Legend: color per sorbent + grey textures per component.
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import os

IL_XLSX    = r'NOAK_Assessment.xlsx'
WAT_XLSX   = r'NOAK_Assessment_Water.xlsx'
OUTPUT_DIR = r'.'
DPI        = 200

# Read FOAK component values from LCCC breakdown spreadsheet
lccc_df = pd.read_excel('LCCCBreakdown.xlsx', header=0, index_col=0)

FOAK_IL_CAPEX   = float(lccc_df.loc['Annualized Capital', 'IL ($/kg CO2)'])
FOAK_IL_ENERGY  = (float(lccc_df.loc['Electricity',    'IL ($/kg CO2)']) +
                   float(lccc_df.loc['Cooling Water',  'IL ($/kg CO2)']) +
                   float(lccc_df.loc['Refrigeration',  'IL ($/kg CO2)']))
FOAK_IL_SORBENT = float(lccc_df.loc['Raw Materials (Sorbent)', 'IL ($/kg CO2)'])
FOAK_IL_FIXED   = (float(lccc_df.loc['LO+M',                    'IL ($/kg CO2)']) +
                   float(lccc_df.loc['Property Tax + Insurance', 'IL ($/kg CO2)']))
FOAK_IL_TOTAL   = float(lccc_df['IL ($/kg CO2)'].sum())

FOAK_WAT_CAPEX  = float(lccc_df.loc['Annualized Capital', 'Water ($/kg CO2)'])
FOAK_WAT_ENERGY = (float(lccc_df.loc['Electricity',    'Water ($/kg CO2)']) +
                   float(lccc_df.loc['Cooling Water',  'Water ($/kg CO2)']) +
                   float(lccc_df.loc['Refrigeration',  'Water ($/kg CO2)']))
FOAK_WAT_FIXED  = (float(lccc_df.loc['LO+M',                    'Water ($/kg CO2)']) +
                   float(lccc_df.loc['Property Tax + Insurance', 'Water ($/kg CO2)']))
FOAK_WAT_TOTAL  = float(lccc_df['Water ($/kg CO2)'].sum())

WATER_BASE = '#2166AC'
IL_BASE    = '#E48B17'

# Component shades — water blues, IL oranges
WAT_SHADES   = ['#2166AC', '#4393C3', '#92C5DE']
IL_SHADES    = ['#8B4513', '#CD853F', '#E48B17', '#FEE090']
COMP_HATCHES = ['', '//', 'xx', '..']
WAT_HATCHES  = ['', '//', '..']   # CAPEX, Energy, Fixed — skips IL Sorbent slot
WAT_COLS = ['NOAK CAPEX', 'NOAK Energy OPEX', 'NOAK Fixed OPEX']
IL_COLS  = ['NOAK CAPEX', 'NOAK Energy OPEX', 'NOAK IL Sorbent', 'NOAK Fixed OPEX']
COMP_NAMES = ['Annualized CAPEX', 'Energy OPEX', 'IL Sorbent', 'Fixed OPEX']

plt.rcParams.update({
    'font.family': 'sans-serif', 'font.size': 13,
    'figure.facecolor': 'white', 'axes.facecolor': 'white',
})

il_df = pd.read_excel(IL_XLSX, header=0)
il_df.columns = il_df.columns.str.strip()
il_df['Scenario'] = il_df['Scenario'].str.strip()
il_df['Il Cost ($/kg)'] = il_df['Il Cost ($/kg)'].astype(float)

wat_df = pd.read_excel(WAT_XLSX, header=0)
wat_df.columns = wat_df.columns.str.strip()
wat_df['Scenario'] = wat_df['Scenario'].str.strip()

scenarios = ['Optimistic', 'Moderate', 'Pessimistic']
il_costs  = [1.0, 3.0, 10.0]
n_s  = len(scenarios)
n_il = len(il_costs)

bar_w       = 0.22
n_per_group = 1 + n_il   # water + 3 IL
group_gap   = 0.35
foak_sep    = 0.6

# FOAK: water at 0, IL at bar_w+0.05
x_foak_wat = 0.0
x_foak_il  = bar_w + 0.05

# NOAK groups start after FOAK
x_start = x_foak_il + foak_sep
group_starts = [x_start + i*(n_per_group*bar_w + group_gap) for i in range(n_s)]

fig, ax = plt.subplots(figsize=(15, 8))

def draw_stack(ax, xpos, colors, hatches, vals, width, total, fs=11):
    bottom = 0.0
    for j, val in enumerate(vals):
        ax.bar(xpos, val, width, bottom=bottom,
               color=colors[j], hatch=hatches[j],
               edgecolor='white', linewidth=0.5, zorder=3)
        bottom += val
    ax.text(xpos, total*1.03, f'${total:,.0f}',
            ha='center', va='bottom', fontsize=fs, fontweight='bold')

# FOAK Water
draw_stack(ax, x_foak_wat,
           WAT_SHADES, WAT_HATCHES,
           [FOAK_WAT_CAPEX, FOAK_WAT_ENERGY, FOAK_WAT_FIXED],
           bar_w*1.2, FOAK_WAT_TOTAL)

# FOAK IL
draw_stack(ax, x_foak_il,
           IL_SHADES, COMP_HATCHES,
           [FOAK_IL_CAPEX, FOAK_IL_ENERGY, FOAK_IL_SORBENT, FOAK_IL_FIXED],
           bar_w*1.2, FOAK_IL_TOTAL)

# NOAK groups — water first then IL scenarios
for si, scen in enumerate(scenarios):
    gs = group_starts[si]

    # Water bar
    wrow = wat_df[wat_df['Scenario']==scen].iloc[0]
    draw_stack(ax, gs,
               WAT_SHADES, WAT_HATCHES,
               [float(wrow[c]) for c in WAT_COLS],
               bar_w, float(wrow['NOAK LCCC']), fs=9)

    # IL bars
    for ii, ilc in enumerate(il_costs):
        irow = il_df[(il_df['Scenario']==scen) &
                     (il_df['Il Cost ($/kg)']==ilc)].iloc[0]
        xpos = gs + (ii+1)*bar_w
        draw_stack(ax, xpos,
                   IL_SHADES, COMP_HATCHES,
                   [float(irow[c]) for c in IL_COLS],
                   bar_w, float(irow['NOAK LCCC']), fs=9)

# Axis
ax.set_yscale('log')
ax.set_ylim(1, FOAK_IL_TOTAL*5)
ax.set_ylabel('LCCC ($/kg CO₂ captured)', fontsize=13)
ax.yaxis.grid(True, which='both', color='#CCCCCC', linestyle='--',
              linewidth=0.5, alpha=0.7, zorder=0)
ax.set_axisbelow(True)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['bottom'].set_visible(False)

# Main x tick labels
foak_center = (x_foak_wat + x_foak_il) / 2
all_x     = [foak_center]
all_lbls  = ['FOAK\nBaseline']
for si, scen in enumerate(scenarios):
    gs = group_starts[si]
    center = gs + (n_per_group-1)*bar_w/2
    all_x.append(center)
    all_lbls.append(f'NOAK\n{scen}')
ax.set_xticks(all_x)
ax.set_xticklabels([''] * len(all_x))
ax.tick_params(axis='x', length=0)

from matplotlib.transforms import blended_transform_factory
trans = blended_transform_factory(ax.transData, ax.transAxes)

# Group labels (lower row)
for xi, lbl in zip(all_x, all_lbls):
    ax.text(xi, -0.14, lbl, ha='center', va='top', fontsize=12,
            transform=trans, color='black')

# Sub-labels (upper row, just below x-axis)
ax.text(x_foak_wat, -0.03, 'Water', ha='center', va='top', fontsize=9, color='black', transform=trans)
ax.text(x_foak_il,  -0.03, 'IL',    ha='center', va='top', fontsize=9, color='black', transform=trans)
for si, scen in enumerate(scenarios):
    gs = group_starts[si]
    ax.text(gs, -0.03, 'W', ha='center', va='top', fontsize=9, color='black', transform=trans)
    for ii, ilc in enumerate(il_costs):
        ax.text(gs+(ii+1)*bar_w, -0.03, f'IL\n${ilc}/kg',
                ha='center', va='top', fontsize=9, color='black', transform=trans)

ax.set_xlim(x_foak_wat - bar_w*0.8,
            group_starts[-1] + n_per_group*bar_w + bar_w*0.5)

# Legend — color per sorbent + grey textures per component
sorbent_handles = [
    mpatches.Patch(facecolor=IL_BASE,    edgecolor='white', label='IL ([bmim][DCA])'),
    mpatches.Patch(facecolor=WATER_BASE, edgecolor='white', label='Water'),
]
comp_handles = [
    mpatches.Patch(facecolor='#888888', hatch='',   edgecolor='#444444', linewidth=0.8, label='Annualized CAPEX'),
    mpatches.Patch(facecolor='#888888', hatch='//', edgecolor='#444444', linewidth=0.8, label='Energy OPEX'),
    mpatches.Patch(facecolor='#888888', hatch='xx', edgecolor='#444444', linewidth=0.8, label='IL Sorbent (IL only)'),
    mpatches.Patch(facecolor='#888888', hatch='..', edgecolor='#444444', linewidth=0.8, label='Fixed OPEX'),
]
ax.legend(handles=sorbent_handles + comp_handles,
          fontsize=11, loc='upper right', frameon=True,
          framealpha=0.95, edgecolor='#cccccc',
          title='Sorbent / Component', title_fontsize=12, ncol=2)

plt.subplots_adjust(bottom=0.22, top=0.90, left=0.07, right=0.98)
fname = os.path.join(OUTPUT_DIR, 'Combined_FOAK_NOAK.png')
plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {fname}')

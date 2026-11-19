"""
generate_SA_figures_updated.py
==============================
Generates one combined figure per SA variable containing:
  - Top row: line plots of raw output values, water and IL overlaid
  - Bottom row: two spider/radar charts (water left, IL right)
    with normalized MCDA scores (5 criteria)

Usage: python generate_SA_figures_updated.py
Requirements: pip install pandas openpyxl matplotlib numpy
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

# ============================================================
# USER SETTINGS
# ============================================================
XLSX_PATH  = r'MCDA_SA.xlsx'
OUTPUT_DIR = r'.'
DPI        = 200

# ============================================================
# COLORS AND STYLE
# ============================================================
WATER_COLOR = '#2166AC'
IL_COLOR    = '#E48B17'
WATER_HATCH = ''
IL_HATCH    = '//'
BASELINE_LS = '--'
NONBASE_LS  = '-'
GRID_COLOR  = '#CCCCCC'

plt.rcParams.update({
    'font.family':      'sans-serif',
    'font.size':        13,
    'axes.titlesize':   13,
    'axes.labelsize':   12,
    'xtick.labelsize':  11.5,
    'ytick.labelsize':  11.5,
    'figure.facecolor': 'white',
    'axes.facecolor':   'white',
    'axes.grid':        True,
    'grid.color':       GRID_COLOR,
    'grid.linewidth':   0.5,
    'grid.alpha':       0.7,
})

# ============================================================
# SA VARIABLE DEFINITIONS
# Sheet2 row indices (0-based, data rows only — skip header rows)
# Water blocks: rows 2-25, IL blocks: rows 28-51
# ============================================================
SA_DEFS = [
    {
        'key':   'absorption_pressure',
        'title': 'Absorption Pressure ',
        'x_label': 'Absorption Pressure (atm)',
        'discrete': False,
        'sheet1_scenarios_water': [
            '3atmAbsorption','6atmAbsorption','Baseline',
            '15atmAbsorption','20atmAbsorption'
        ],
        'x_tick_labels': ['3','6','10\n(BL)','15','20'],
        'baseline_idx': 2,
        'sheet2_water_start': 2,
        'sheet2_water_end':   7,
        'sheet2_il_start':    28,
        'sheet2_il_end':      33,
    },
    {
        'key':   'velocity',
        'title': 'Superficial Gas Velocity',
        'x_label': 'Superficial Gas Velocity (m/s)',
        'discrete': False,
        'sheet1_scenarios_water': [
            '0.005mpsVelocity','0.01mpsVelocity','0.025mpsVelocity',
            'Baseline','0.1mpsVelocity','0.2mpsVelocity'
        ],
        'x_tick_labels': ['0.005','0.01','0.025','0.05\n(BL)','0.1','0.2'],
        'baseline_idx': 3,
        'sheet2_water_start': 8,
        'sheet2_water_end':   14,
        'sheet2_il_start':    34,
        'sheet2_il_end':      40,
    },
    {
        'key':   'temperature',
        'title': 'Absorber Temperature',
        'x_label': 'Absorber Temperature (°C)',
        'discrete': False,
        'sheet1_scenarios_water': ['21CTemp','Baseline','30CTemp'],
        'x_tick_labels': ['21','25\n(BL)','30'],
        'baseline_idx': 1,
        'sheet2_water_start': 15,
        'sheet2_water_end':   18,
        'sheet2_il_start':    41,
        'sheet2_il_end':      44,
    },
    {
        'key':   'desorption_pressure',
        'title': 'Vacuum Desorption Pressure',
        'x_label': 'Desorption Pressure (atm)',
        'discrete': False,
        'sheet1_scenarios_water': [
            '0.05atmDesorption','Baseline','0.2atmDesorption'
        ],
        'x_tick_labels': ['0.05','0.1\n(BL)','0.2'],
        'baseline_idx': 1,
        'sheet2_water_start': 19,
        'sheet2_water_end':   22,
        'sheet2_il_start':    45,
        'sheet2_il_end':      48,
    },
    {
        'key':   'purification_temp',
        'title': 'CO\u2082 Purification Temperature',
        'x_label': 'Purification Temperature (°C)',
        'discrete': False,
        'sheet1_scenarios_water': [
            '-30CPurification','Baseline','-90CPurification'
        ],
        'x_tick_labels': ['-30','-60\n(BL)','-90'],
        'baseline_idx': 1,
        'sheet2_water_start': 23,
        'sheet2_water_end':   26,
        'sheet2_il_start':    49,
        'sheet2_il_end':      52,
    },
]

# Raw metric columns from Sheet1 — 5 metrics
RAW_METRICS = [
    ('Capture Efficiency (%)',       'Capture Efficiency (%)'),
    ('CO2Product (kg/h)',             'Capture Rate (kg/hr)'),
    ('CO2ProductPurity (wt%)',        'Product Purity (wt%)'),
    ('Net Electricity (MWh/kg CO2)', 'Net Electricity (MWh/kg CO\u2082)'),
    ('Sorbent Makeup (kg/hr)',        'Sorbent Makeup (kg/hr)'),
]

# Spider axis labels — 5 criteria
SPIDER_LABELS = [
    'Capture\nEfficiency',
    'Capture\nRate',
    'Product\nPurity',
    'Energy',
    'Makeup\nRate',
]

# Sheet2 normalized score columns (0-based)
# CE norm=2, CR norm=4, PP norm=6, EC norm=8, MR norm=10
SPIDER_SCORE_COLS = [2, 4, 6, 8, 10]

# ============================================================
# LOAD DATA
# ============================================================
print(f'Reading {XLSX_PATH} ...')
df1_raw = pd.read_excel(XLSX_PATH, sheet_name='Sheet1', header=0)
df2_raw = pd.read_excel(XLSX_PATH, sheet_name='Sheet2', header=None)

df1_water = df1_raw[df1_raw['Sorbent'] == 'Water'].copy().reset_index(drop=True)
df1_il    = df1_raw[df1_raw['Sorbent'] == 'IL'].copy().reset_index(drop=True)

# Add Baseline row mapping
def get_sheet1_rows(df, scenario_list):
    rows = []
    for sc in scenario_list:
        if sc == 'Baseline':
            match = df[df['Scenario'] == 'Baseline']
        else:
            match = df[df['Scenario'] == sc]
        if len(match) == 0:
            raise ValueError(f"Scenario '{sc}' not found.")
        rows.append(match.iloc[0])
    return pd.DataFrame(rows).reset_index(drop=True)

def get_sheet2_block(start, end):
    return df2_raw.iloc[start:end].reset_index(drop=True)

def build_spider_dict(block, x_labels):
    d = {}
    for i, lbl in enumerate(x_labels):
        scores = [float(block.iloc[i, c]) for c in SPIDER_SCORE_COLS]
        d[lbl] = scores
    return d

# ============================================================
# SPIDER PLOT HELPER
# ============================================================
def draw_spider(ax, data_dict, colors, linestyles, linewidths, title):
    N = len(SPIDER_LABELS)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles += angles[:1]

    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_xticks(angles[:-1])
    ax.set_xticklabels([])
    ax.set_ylim(0, 1)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(['0.25','0.50','0.75','1.00'],
                       color='grey')
    ax.tick_params(pad=9)
    ax.grid(color='grey', linestyle='--', linewidth=0.4, alpha=0.5)
    ax.spines['polar'].set_visible(False)
    ax.fill(angles, [1]*len(angles), alpha=0.04, color='grey')

    label_radius = 1.05
    for angle, label in zip(angles[:-1], SPIDER_LABELS):
        ha = ('center' if np.isclose(angle, 0) or np.isclose(angle, np.pi)
            else 'left' if 0 < angle < np.pi else 'right')
        va = ('center' if np.isclose(angle, np.pi/2) or
            np.isclose(angle, 3*np.pi/2)
            else 'bottom' if angle < np.pi/2 or angle > 3*np.pi/2
            else 'top')
        # Push top label further out to avoid title overlap
        radius = 1.20 if np.isclose(angle, np.pi/2) else label_radius
        ax.text(angle, radius, label,
                ha=ha, va=va, multialignment='center')

    handles = []
    for i, (lbl, scores) in enumerate(data_dict.items()):
        vals = list(scores) + [scores[0]]
        ax.plot(angles, vals, color=colors[i], linewidth=linewidths[i],
                linestyle=linestyles[i], alpha=0.9, zorder=3)
        ax.fill(angles, vals, color=colors[i], alpha=0.10)
        handles.append(mpatches.Patch(
            color=colors[i],
            label=lbl.replace('\n', ' '),
            alpha=0.85))

    ax.set_title(title, fontweight='bold', pad=55, fontsize=12)
    return handles

# ============================================================
# MAIN FIGURE GENERATOR
# ============================================================
def make_figure(sa):
    key      = sa['key']
    n_pts    = len(sa['x_tick_labels'])
    bl_idx   = sa['baseline_idx']
    x_labels = sa['x_tick_labels']
    x        = np.arange(n_pts)

    w_rows  = get_sheet1_rows(df1_water, sa['sheet1_scenarios_water'])
    il_rows = get_sheet1_rows(df1_il,    sa['sheet1_scenarios_water'])

    w2 = get_sheet2_block(sa['sheet2_water_start'], sa['sheet2_water_end'])
    i2 = get_sheet2_block(sa['sheet2_il_start'],    sa['sheet2_il_end'])
    w_spider = build_spider_dict(w2, x_labels)
    i_spider = build_spider_dict(i2, x_labels)

    # Spider colors
    import matplotlib.colors as mcolors
    cmap_w   = plt.cm.Blues(np.linspace(0.40, 0.90, n_pts))
    w_colors = [cmap_w[j] for j in range(n_pts)]
    base_rgb = mcolors.to_rgb(IL_COLOR)
    i_colors = [
        tuple(1-(1-c)*f for c in base_rgb)
        for f in np.linspace(0.30, 1.0, n_pts)
    ]
    w_ls = [BASELINE_LS if j == bl_idx else NONBASE_LS for j in range(n_pts)]
    i_ls = w_ls[:]
    w_lw = [2.5 if j == bl_idx else 1.6 for j in range(n_pts)]
    i_lw = w_lw[:]

    n_metrics = len(RAW_METRICS)
    fig = plt.figure(figsize=(n_metrics * 3.4, 12))
    gs  = gridspec.GridSpec(
        2, n_metrics,
        figure=fig,
        height_ratios=[1.1, 0.75],
        hspace=0.6,
        wspace=0.55,
        top=0.93,
        bottom=0.12,
        left=0.06,
        right=0.98,
    )

    bar_width = 0.35

    for m_idx, (col, ylabel) in enumerate(RAW_METRICS):
        ax = fig.add_subplot(gs[0, m_idx])

        w_vals  = w_rows[col].values.astype(float)
        il_vals = il_rows[col].values.astype(float)

        # Handle missing -30C water purification
        if key == 'purification_temp':
            if col == 'Capture Efficiency (%)' and w_vals[0] == 0:
                w_vals[0] = np.nan

        if sa['discrete']:
            ax.bar(x - bar_width/2, w_vals,  bar_width,
                   color=WATER_COLOR, hatch=WATER_HATCH,
                   edgecolor='white', alpha=0.85, label='Water')
            ax.bar(x + bar_width/2, il_vals, bar_width,
                   color=IL_COLOR,   hatch=IL_HATCH,
                   edgecolor='white', alpha=0.85,
                   label='IL ([bmim][DCA])')
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels)
        else:
            ax.plot(x, w_vals,  color=WATER_COLOR, marker='o',
                    linewidth=1.8, markersize=5,
                    label='Water', zorder=3)
            ax.plot(x, il_vals, color=IL_COLOR,    marker='s',
                    linewidth=1.8, markersize=5, linestyle='--',
                    label='IL ([bmim][DCA])', zorder=3)
            ax.axvline(bl_idx, color='grey', linestyle=':',
                       linewidth=1.0, alpha=0.6, zorder=1)
            ax.set_xticks(x)
            ax.set_xticklabels(x_labels)

        ax.set_ylabel(ylabel)
        ax.set_xlabel(sa['x_label'])
        ax.set_title(ylabel, fontweight='bold')

        if m_idx == 0:
            ax.legend(fontsize=10, loc='best', framealpha=0.9)

    # Spider charts — use a sub-GridSpec so both halves are exactly equal width
    from matplotlib.gridspec import GridSpecFromSubplotSpec
    gs_bottom = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, :], wspace=0.4)
    ax_sw = fig.add_subplot(gs_bottom[0], polar=True)
    ax_si = fig.add_subplot(gs_bottom[1], polar=True)

    h_w  = draw_spider(ax_sw, w_spider, w_colors, w_ls, w_lw,
                       'Water Sorbent')
    h_il = draw_spider(ax_si, i_spider, i_colors, i_ls, i_lw,
                       'IL Sorbent ([bmim][DCA])')

    ax_sw.legend(handles=h_w,  loc='upper center',
                 bbox_to_anchor=(0.5, -0.15), ncol=2,
                 fontsize=10, frameon=True, framealpha=0.9,
                 edgecolor='#cccccc',
                 title='Setpoint (dashed = baseline)',
                 title_fontsize=10)
    ax_si.legend(handles=h_il, loc='upper center',
                 bbox_to_anchor=(0.5, -0.15), ncol=2,
                 fontsize=10, frameon=True, framealpha=0.9,
                 edgecolor='#cccccc',
                 title='Setpoint (dashed = baseline)',
                 title_fontsize=10)


    fname = os.path.join(OUTPUT_DIR, f'SA_{key}.png')
    plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {fname}')

# ============================================================
# TORNADO PLOT
# ============================================================
def _clean_tick(lbl):
    """Strip newline and (BL) suffix from x_tick_labels entries."""
    return lbl.replace('\n(BL)', '').replace('\n', '').strip()


def make_tornado():
    # Build y-axis labels: "Parameter Name\n(low, Baseline, high)"
    sa_labels = []
    for sa in SA_DEFS:
        ticks  = sa['x_tick_labels']
        bl_idx = sa['baseline_idx']
        low    = _clean_tick(ticks[0])
        base   = _clean_tick(ticks[bl_idx])
        high   = _clean_tick(ticks[-1])
        sa_labels.append(f"{sa['title']}\n({low}, $\\mathbf{{{base}}}$, {high})")
    sa_scenarios = [
        ['3atmAbsorption','6atmAbsorption','Baseline','15atmAbsorption','20atmAbsorption'],
        ['0.005mpsVelocity','0.01mpsVelocity','0.025mpsVelocity','Baseline','0.1mpsVelocity','0.2mpsVelocity'],
        ['21CTemp','Baseline','30CTemp'],
        ['0.05atmDesorption','Baseline','0.2atmDesorption'],
        ['-30CPurification','Baseline','-90CPurification'],
    ]

    metrics = [
        ('Net Electricity (MWh/kg CO2)', 'Specific Energy\n(MWh/kg CO\u2082)'),
        ('Capture Efficiency (%)',        'Capture Efficiency (%)'),
    ]

    n_sa     = len(sa_labels)
    y_pos    = np.arange(n_sa)
    height   = 0.35

    fig, axes = plt.subplots(1, 2, figsize=(14, 7), sharey=True)

    for ax, (col, xlabel) in zip(axes, metrics):
        for i, (scenarios, _) in enumerate(zip(sa_scenarios, sa_labels)):
            w_sc = scenarios
            w_vals  = get_sheet1_rows(df1_water, w_sc)[col].values.astype(float)
            w_base  = get_sheet1_rows(df1_water, ['Baseline'])[col].values[0]
            w_min   = w_vals.min() - w_base
            w_max   = w_vals.max() - w_base

            il_vals = get_sheet1_rows(df1_il, scenarios)[col].values.astype(float)
            il_base = get_sheet1_rows(df1_il, ['Baseline'])[col].values[0]
            il_min  = il_vals.min() - il_base
            il_max  = il_vals.max() - il_base

            ax.barh(y_pos[i] + height/2, w_max - w_min,
                    left=w_min, height=height,
                    color=WATER_COLOR, alpha=0.85,
                    label='Water' if i == 0 else '')
            ax.barh(y_pos[i] - height/2, il_max - il_min,
                    left=il_min, height=height,
                    color=IL_COLOR, hatch='//', alpha=0.85,
                    label='IL ([bmim][DCA])' if i == 0 else '')

        ax.axvline(0, color='black', linewidth=1.0,
                   linestyle='--', alpha=0.6)
        xabs = max(abs(ax.get_xlim()[0]), abs(ax.get_xlim()[1]))
        ax.set_xlim(-xabs, xabs)
        ax.set_xlabel(f'\u0394 {xlabel}\n(relative to baseline)')
        ax.set_title(xlabel, fontweight='bold')

    axes[0].set_yticks(y_pos)
    axes[0].set_yticklabels(sa_labels)
    axes[0].set_ylim(n_sa - 0.5, -0.5)
    axes[1].set_ylim(n_sa - 0.5, -0.5)
    axes[0].legend(fontsize=11, loc='lower right', framealpha=0.9)

    plt.tight_layout()
    fname = os.path.join(OUTPUT_DIR, 'SA_tornado.png')
    plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f'  Saved: {fname}')

# ============================================================
# RUN
# ============================================================
if __name__ == '__main__':
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    for sa_def in SA_DEFS:
        print(f'Generating: {sa_def["title"]} ...')
        try:
            make_figure(sa_def)
        except Exception as e:
            print(f'  ERROR: {e}')
    print('Generating: Tornado plot ...')
    try:
        make_tornado()
    except Exception as e:
        print(f'  ERROR: {e}')
    print('\nDone.')

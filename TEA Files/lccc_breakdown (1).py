"""
lccc_breakdown.py
=================
Stacked bar chart of LCCC breakdown with broken y-axis.
Bottom panel: 0 to BREAK_LOW (shows all components clearly)
Top panel: BREAK_HIGH to max (shows IL sorbent cost)

Usage: python lccc_breakdown.py
Requirements: pip install pandas openpyxl matplotlib numpy
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
XLSX_PATH  = r'C:\Users\anbarron\Box\Research\Projects\DAC\Modeling\TEA\LCCCBreakdown.xlsx'
OUTPUT_DIR = r'C:\Users\anbarron\Box\Research\Projects\DAC\Modeling\TEA'
DPI        = 200

BREAK_LOW  = 225    # top of lower panel ($/kg CO2); must exceed water total (202)
BREAK_HIGH = 3800   # bottom of upper panel ($/kg CO2)

# ============================================================
# COLORS
# ============================================================
WATER_BASE = '#2166AC'
IL_BASE    = '#E48B17'
HATCHES    = ['', '//', '\\\\', 'xx', '..', '++', '--']

def make_shades(base_hex, n):
    import matplotlib.colors as mc
    rgb = mc.to_rgb(base_hex)
    return [tuple(1-(1-c)*f for c in rgb)
            for f in np.linspace(0.40, 0.92, n)]

# ============================================================
# LOAD DATA
# ============================================================
df           = pd.read_excel(XLSX_PATH, header=0, index_col=0)
df.columns   = ['Water', 'IL']
components   = df.index.tolist()
n_comp       = len(components)
water_vals   = df['Water'].values.astype(float)
il_vals      = df['IL'].values.astype(float)
water_shades = make_shades(WATER_BASE, n_comp)
il_shades    = make_shades(IL_BASE,    n_comp)

# Reorder so Raw Materials plots last (top of stack)
desired_order = [
    'Annualized Capital',
    'LO+M',
    'Cooling Water',
    'Electricity',
    'Property Tax + Insurance',
    'Refrigeration',
    'Raw Materials (Sorbent)',
]
df = df.reindex(desired_order)
components   = df.index.tolist()
water_vals   = df['Water'].values.astype(float)
il_vals      = df['IL'].values.astype(float)
water_shades = make_shades(WATER_BASE, n_comp)
il_shades    = make_shades(IL_BASE,    n_comp)

# ============================================================
# FIGURE — left: broken-axis LCCC; right: 100% contribution
# ============================================================
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(13, 8))
fig.patch.set_facecolor('white')

gs = gridspec.GridSpec(2, 2,
                       height_ratios=[1, 2.5],
                       hspace=0.06, wspace=0.45,
                       figure=fig)
ax_top = fig.add_subplot(gs[0, 0])
ax_bot = fig.add_subplot(gs[1, 0], sharex=ax_top)
ax_pct = fig.add_subplot(gs[:, 1])

for ax in (ax_top, ax_bot, ax_pct):
    ax.set_facecolor('white')

bar_width = 0.35
x = np.array([0, 1])

# ---- Draw stacked bars on BOTH axes -----------------------
water_bottom = 0.0
il_bottom    = 0.0
legend_handles = []

for i, comp in enumerate(components):
    hatch = HATCHES[i % len(HATCHES)]

    for ax in (ax_top, ax_bot):
        ax.bar(x[0], water_vals[i], bar_width,
               bottom=water_bottom,
               color=water_shades[i],
               hatch=hatch,
               edgecolor='white',
               linewidth=0.5,
               zorder=3)
        ax.bar(x[1], il_vals[i], bar_width,
               bottom=il_bottom,
               color=il_shades[i],
               hatch=hatch,
               edgecolor='white',
               linewidth=0.5,
               zorder=3)

    water_bottom += water_vals[i]
    il_bottom    += il_vals[i]

    legend_handles.append(mpatches.Patch(
        facecolor='#888888',
        hatch=hatch,
        edgecolor='#444444',
        linewidth=0.8,
        label=comp
    ))

# ---- Set axis limits for broken effect --------------------
ax_bot.set_ylim(0, BREAK_LOW)
ax_top.set_ylim(BREAK_HIGH, il_bottom * 1.10)

# ---- Total annotations ------------------------------------
ax_top.text(x[1], il_bottom * 1.02,
            f'${il_bottom:.0f}/kg CO₂',
            ha='center', va='bottom', fontsize=12,
            fontweight='bold', color='black')
ax_bot.text(x[0], water_bottom + BREAK_LOW * 0.02,
            f'${water_bottom:.0f}/kg CO₂',
            ha='center', va='bottom', fontsize=12,
            fontweight='bold', color='black')

# ---- Hide inner spines and add break marks ----------------
ax_top.spines['bottom'].set_visible(False)
ax_bot.spines['top'].set_visible(False)
ax_top.tick_params(bottom=False, labeltop=False)
ax_top.xaxis.tick_top()
ax_top.set_xticks([])

# Diagonal break lines
d = 0.012
kwargs = dict(transform=ax_top.transAxes, color='k',
              linewidth=1.2, clip_on=False)
ax_top.plot((-d, +d), (-d, +d), **kwargs)
ax_top.plot((1-d, 1+d), (-d, +d), **kwargs)

kwargs.update(transform=ax_bot.transAxes)
ax_bot.plot((-d, +d), (1-d, 1+d), **kwargs)
ax_bot.plot((1-d, 1+d), (1-d, 1+d), **kwargs)

# ---- Grid -------------------------------------------------
for ax in (ax_top, ax_bot):
    ax.yaxis.grid(True, color='#CCCCCC', linestyle='--',
                  linewidth=0.5, alpha=0.7, zorder=0)
    ax.set_axisbelow(True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

# ---- Labels (left panel) ----------------------------------
ax_bot.set_xticks(x)
ax_bot.set_xticklabels(['Water Sorbent', 'IL Sorbent\n([bmim][DCA])'],
                        fontsize=13)
plt.setp(ax_top.get_xticklabels(), visible=False)
fig.text(0.02, 0.5, 'LCCC ($/kg CO₂ captured)',
         va='center', rotation='vertical', fontsize=13)

# ---- 100% contribution bar chart (right panel) ------------
water_total = water_bottom
il_total    = il_bottom

w_bot_pct = 0.0
i_bot_pct = 0.0
for i, comp in enumerate(components):
    hatch  = HATCHES[i % len(HATCHES)]
    wv_pct = water_vals[i] / water_total * 100
    iv_pct = il_vals[i]    / il_total    * 100
    ax_pct.bar(x[0], wv_pct, bar_width,
               bottom=w_bot_pct,
               color=water_shades[i], hatch=hatch,
               edgecolor='white', linewidth=0.5, zorder=3)
    ax_pct.bar(x[1], iv_pct, bar_width,
               bottom=i_bot_pct,
               color=il_shades[i], hatch=hatch,
               edgecolor='white', linewidth=0.5, zorder=3)
    w_bot_pct += wv_pct
    i_bot_pct += iv_pct

ax_pct.set_ylim(0, 100)
ax_pct.set_xticks(x)
ax_pct.set_xticklabels(['Water Sorbent', 'IL Sorbent\n([bmim][DCA])'],
                        fontsize=13)
ax_pct.set_ylabel('Contribution to LCCC (%)', fontsize=13)
ax_pct.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
ax_pct.yaxis.grid(True, color='#CCCCCC', linestyle='--',
                  linewidth=0.5, alpha=0.7, zorder=0)
ax_pct.set_axisbelow(True)
ax_pct.spines['top'].set_visible(False)
ax_pct.spines['right'].set_visible(False)

# ---- Shared legend ----------------------------------------
water_key = mpatches.Patch(facecolor=WATER_BASE, edgecolor='white',
                            label='Water (blue tones)')
il_key    = mpatches.Patch(facecolor=IL_BASE,    edgecolor='white',
                            label='IL (orange tones)')
fig.legend(
    handles=legend_handles + [water_key, il_key],
    title='Cost Component',
    title_fontsize=12,
    fontsize=11,
    loc='lower center',
    bbox_to_anchor=(0.5, -0.12),
    ncol=3,
    frameon=True,
    framealpha=0.95,
    edgecolor='#cccccc',
    handlelength=3.5,
    handleheight=2.0,
)

# ---- Panel labels (figure coords so (a)/(b) share same y) -
plt.tight_layout()
pos_top = ax_top.get_position()
pos_pct = ax_pct.get_position()
label_y = pos_top.y1 + 0.02
fig.text(pos_top.x0 - 0.01, label_y, '(a)',
         fontsize=14, fontweight='bold', va='bottom', ha='right',
         transform=fig.transFigure)
fig.text(pos_pct.x0 - 0.01, label_y, '(b)',
         fontsize=14, fontweight='bold', va='bottom', ha='right',
         transform=fig.transFigure)

# ---- Save -------------------------------------------------
fname = os.path.join(OUTPUT_DIR, 'LCCC_breakdown.png')
plt.savefig(fname, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved: {fname}')

# ============================================================
# FIGURE 2 — left: log-scale LCCC (no break); right: 100% contribution
# ============================================================
fig2, (ax_log, ax_pct2) = plt.subplots(1, 2, figsize=(13, 7))
fig2.patch.set_facecolor('white')

for ax in (ax_log, ax_pct2):
    ax.set_facecolor('white')

# ---- Log-scale stacked bars -------------------------------
water_bottom2 = 1e-9   # small nonzero base for log stacking
il_bottom2    = 1e-9
water_bottoms2 = []
il_bottoms2    = []

for val in water_vals:
    water_bottoms2.append(water_bottom2)
    water_bottom2 += val
for val in il_vals:
    il_bottoms2.append(il_bottom2)
    il_bottom2 += val

for i, comp in enumerate(components):
    hatch = HATCHES[i % len(HATCHES)]
    ax_log.bar(x[0], water_vals[i], bar_width,
               bottom=water_bottoms2[i],
               color=water_shades[i], hatch=hatch,
               edgecolor='white', linewidth=0.5, zorder=3)
    ax_log.bar(x[1], il_vals[i], bar_width,
               bottom=il_bottoms2[i],
               color=il_shades[i], hatch=hatch,
               edgecolor='white', linewidth=0.5, zorder=3)

ax_log.set_yscale('log')
ax_log.set_ylim(bottom=min(v for v in np.concatenate([water_vals, il_vals]) if v > 0) * 0.5)
ax_log.set_xticks(x)
ax_log.set_xticklabels(['Water Sorbent', 'IL Sorbent\n([bmim][DCA])'], fontsize=13)
ax_log.set_ylabel('LCCC ($/kg CO\u2082 captured)', fontsize=13)
ax_log.yaxis.grid(True, color='#CCCCCC', linestyle='--',
                  linewidth=0.5, alpha=0.7, zorder=0)
ax_log.set_axisbelow(True)
ax_log.spines['top'].set_visible(False)
ax_log.spines['right'].set_visible(False)

# Total annotations
ax_log.text(x[0], water_bottom2 * 1.15,
            f'${water_bottom2:.0f}/kg CO\u2082',
            ha='center', va='bottom', fontsize=12,
            fontweight='bold', color='black')
ax_log.text(x[1], il_bottom2 * 1.15,
            f'${il_bottom2:.0f}/kg CO\u2082',
            ha='center', va='bottom', fontsize=12,
            fontweight='bold', color='black')

# ---- 100% contribution bar chart (right panel, fig2) ------
water_total2 = sum(water_vals)
il_total2    = sum(il_vals)
w_bot_pct2 = 0.0
i_bot_pct2 = 0.0

for i, comp in enumerate(components):
    hatch  = HATCHES[i % len(HATCHES)]
    wv_pct = water_vals[i] / water_total2 * 100
    iv_pct = il_vals[i]    / il_total2    * 100
    ax_pct2.bar(x[0], wv_pct, bar_width,
                bottom=w_bot_pct2,
                color=water_shades[i], hatch=hatch,
                edgecolor='white', linewidth=0.5, zorder=3)
    ax_pct2.bar(x[1], iv_pct, bar_width,
                bottom=i_bot_pct2,
                color=il_shades[i], hatch=hatch,
                edgecolor='white', linewidth=0.5, zorder=3)
    w_bot_pct2 += wv_pct
    i_bot_pct2 += iv_pct

ax_pct2.set_ylim(0, 100)
ax_pct2.set_xticks(x)
ax_pct2.set_xticklabels(['Water Sorbent', 'IL Sorbent\n([bmim][DCA])'], fontsize=13)
ax_pct2.set_ylabel('Contribution to LCCC (%)', fontsize=13)
ax_pct2.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.0f}%'))
ax_pct2.yaxis.grid(True, color='#CCCCCC', linestyle='--',
                   linewidth=0.5, alpha=0.7, zorder=0)
ax_pct2.set_axisbelow(True)
ax_pct2.spines['top'].set_visible(False)
ax_pct2.spines['right'].set_visible(False)

# ---- Shared legend ----------------------------------------
fig2.legend(
    handles=legend_handles + [water_key, il_key],
    title='Cost Component',
    title_fontsize=12,
    fontsize=11,
    loc='lower center',
    bbox_to_anchor=(0.5, -0.22),   # lowered so it clears the x tick labels
    ncol=3,
    frameon=True,
    framealpha=0.95,
    edgecolor='#cccccc',
    handlelength=3.5,
    handleheight=2.0,
)

# ---- Panel labels -----------------------------------------
plt.tight_layout()
pos_log  = ax_log.get_position()
pos_pct2 = ax_pct2.get_position()
label_y2 = max(pos_log.y1, pos_pct2.y1) + 0.02
fig2.text(pos_log.x0 - 0.01,  label_y2, '(a)',
          fontsize=14, fontweight='bold', va='bottom', ha='right',
          transform=fig2.transFigure)
fig2.text(pos_pct2.x0 - 0.01, label_y2, '(b)',
          fontsize=14, fontweight='bold', va='bottom', ha='right',
          transform=fig2.transFigure)

# ---- Save -------------------------------------------------
fname2 = os.path.join(OUTPUT_DIR, 'LCCC_breakdown_logscale.png')
fig2.savefig(fname2, dpi=DPI, bbox_inches='tight', facecolor='white')
plt.close(fig2)
print(f'Saved: {fname2}')

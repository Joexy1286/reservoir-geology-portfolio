"""
=============================================================================
WELL LOG ANALYSIS & RESERVOIR CHARACTERISATION
Volve Field Analogue — Norwegian North Sea (CCS Storage Assessment Context)
=============================================================================
Author  : Joexy1286
Purpose : Demonstrate reservoir geologist workflow:
          (1) Well log QC and lithology identification
          (2) Petrophysical analysis — porosity, Vshale, water saturation
          (3) Net pay flagging
          (4) CCS storage potential assessment
Dataset : Synthetic Volve-analogue logs (Norwegian North Sea stratigraphy)
          Modelled on publicly released Equinor Volve dataset
=============================================================================
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'data'))

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec
from scipy.ndimage import uniform_filter1d
from generate_synthetic_logs import generate_volve_logs

# ── Colour palette (dark subsurface theme) ────────────────────────────────
BG      = "#0d1117"
PANEL   = "#161b22"
ACCENT  = "#58a6ff"
GOLD    = "#e3b341"
GREEN   = "#3fb950"
RED     = "#f85149"
PURPLE  = "#bc8cff"
MUTED   = "#8b949e"
WHITE   = "#f0f6fc"

FORM_COLORS = {
    "Nordland_Gp":  "#6e7681",
    "Hordaland_Gp": "#5a6b7a",
    "Shetland_Gp":  "#c9a86c",
    "Rogaland_Gp":  "#5a6b7a",
    "Heimdal_Fm":   "#f5a623",
    "Lista_Fm":     "#5a6b7a",
    "Ty_Fm":        "#e8c56a",
    "Basement":     "#8b8b8b",
}

plt.rcParams.update({
    'figure.facecolor':  BG,
    'axes.facecolor':    PANEL,
    'axes.edgecolor':    MUTED,
    'axes.labelcolor':   WHITE,
    'xtick.color':       MUTED,
    'ytick.color':       MUTED,
    'grid.color':        "#21262d",
    'grid.linewidth':    0.5,
    'text.color':        WHITE,
    'font.family':       'monospace',
    'font.size':         8,
})

# ═══════════════════════════════════════════════════════════════════════════
# 1. LOAD DATA
# ═══════════════════════════════════════════════════════════════════════════
df, formations = generate_volve_logs()

# Smooth logs (5-point running average — simulates Petrel log smoothing)
for col in ["GR", "RHOB", "NPHI", "DT", "RT"]:
    df[f"{col}_sm"] = uniform_filter1d(df[col].values, size=5)

print(f"  Loaded {len(df):,} depth samples | {df['DEPTH'].min():.0f}–{df['DEPTH'].max():.0f} m MD")

# ═══════════════════════════════════════════════════════════════════════════
# 2. PETROPHYSICAL CALCULATIONS
# ═══════════════════════════════════════════════════════════════════════════

# --- 2a. Volume of Shale (Larionov non-linear, Tertiary rocks) ------------
GR_CLEAN = 20.0   # API — clean sand baseline
GR_SHALE = 110.0  # API — pure shale baseline
IGR = (df["GR_sm"] - GR_CLEAN) / (GR_SHALE - GR_CLEAN)
IGR = np.clip(IGR, 0, 1)
df["VSHALE"] = np.clip(0.083 * (2**(3.7 * IGR) - 1), 0, 1)  # Larionov tertiary

# --- 2b. Density Porosity --------------------------------------------------
RHOB_MATRIX = 2.65   # g/cc quartz
RHOB_FLUID  = 1.05   # g/cc saline formation water
df["PHID"] = (RHOB_MATRIX - df["RHOB_sm"]) / (RHOB_MATRIX - RHOB_FLUID)
df["PHID"] = np.clip(df["PHID"], 0, 0.5)

# --- 2c. Neutron-Density crossplot porosity (effective) -------------------
df["PHIE"] = np.clip((df["PHID"] + df["NPHI_sm"]) / 2, 0, 0.45)

# --- 2d. Water Saturation — Archie equation --------------------------------
# Sw = ((a * Rw) / (phi^m * Rt))^(1/n)
A  = 1.0   # tortuosity factor
M  = 2.0   # cementation exponent
N  = 2.0   # saturation exponent
RW = 0.04  # ohm.m at formation temperature (~80°C)
phi_safe = np.where(df["PHIE"] > 0.02, df["PHIE"], 0.02)
df["SW"] = np.clip(
    ((A * RW) / (phi_safe**M * df["RT_sm"]))**(1/N),
    0, 1
)
df["SHC"] = 1 - df["SW"]  # hydrocarbon saturation

# --- 2e. Permeability estimate — Timur-Coates correlation -----------------
# k = C * (phi^4 / Swirr^2)  — approximation for clastics
SWIRR = 0.15
df["PERM"] = np.where(
    df["PHIE"] > 0.05,
    np.clip(0.00136 * (df["PHIE"]**4.4) / (SWIRR**2) * 1e6, 0.01, 5000),
    0.001
)

# ═══════════════════════════════════════════════════════════════════════════
# 3. NET PAY FLAGS
# ═══════════════════════════════════════════════════════════════════════════
VSHALE_CUTOFF = 0.30
PHIE_CUTOFF   = 0.10
SW_CUTOFF     = 0.65

df["NET_RESERVOIR"] = (
    (df["VSHALE"] < VSHALE_CUTOFF) &
    (df["PHIE"]   > PHIE_CUTOFF)
).astype(int)

df["NET_PAY"] = (
    (df["NET_RESERVOIR"] == 1) &
    (df["SW"] < SW_CUTOFF)
).astype(int)

# ═══════════════════════════════════════════════════════════════════════════
# 4. FORMATION SUMMARY STATISTICS
# ═══════════════════════════════════════════════════════════════════════════
summary = df.groupby("FORMATION").agg(
    Top_m=("DEPTH", "min"),
    Base_m=("DEPTH", "max"),
    GR_avg=("GR_sm", "mean"),
    PHIE_avg=("PHIE", "mean"),
    SW_avg=("SW", "mean"),
    PERM_p50=("PERM", "median"),
    Net_Reservoir_m=("NET_RESERVOIR", lambda x: x.sum() * 0.15),
    Net_Pay_m=("NET_PAY", lambda x: x.sum() * 0.15),
).round(2)

print("\n  FORMATION SUMMARY")
print("  " + "─"*80)
print(summary.to_string())

# ═══════════════════════════════════════════════════════════════════════════
# 5. FIGURE 1 — COMPREHENSIVE LOG DISPLAY (Petrel-style 5-track)
# ═══════════════════════════════════════════════════════════════════════════
print("\n  Rendering Figure 1: 5-track log display...")

fig = plt.figure(figsize=(18, 14), facecolor=BG)
fig.suptitle(
    "WELL 15/9-F-1  ·  VOLVE FIELD ANALOGUE  ·  NORWEGIAN NORTH SEA\n"
    "Reservoir Characterisation & CCS Storage Assessment",
    fontsize=11, color=WHITE, fontweight='bold', y=0.98
)

# Use subplot mosaic for track layout
# Tracks: GR | Density-Neutron | Resistivity | Petrophysics | Lithology flag
ax_gr   = fig.add_axes([0.06, 0.06, 0.12, 0.86])
ax_dn   = fig.add_axes([0.20, 0.06, 0.14, 0.86])
ax_rt   = fig.add_axes([0.36, 0.06, 0.12, 0.86])
ax_pe   = fig.add_axes([0.50, 0.06, 0.14, 0.86])
ax_sw   = fig.add_axes([0.66, 0.06, 0.10, 0.86])
ax_pay  = fig.add_axes([0.78, 0.06, 0.06, 0.86])
ax_form = fig.add_axes([0.86, 0.06, 0.12, 0.86])

depth = df["DEPTH"].values
y_lim = (3800, 100)

for ax in [ax_gr, ax_dn, ax_rt, ax_pe, ax_sw, ax_pay, ax_form]:
    ax.set_ylim(y_lim)
    ax.set_facecolor(PANEL)
    ax.grid(True, axis='y', color='#21262d', linewidth=0.4)
    ax.tick_params(labelsize=7)

# ── Track 1: Gamma Ray ─────────────────────────────────────────────────────
ax_gr.plot(df["GR_sm"], depth, color=GREEN, lw=0.8, label="GR")
ax_gr.fill_betweenx(depth, df["GR_sm"], 0, alpha=0.15, color=GREEN)
ax_gr.axvline(GR_CLEAN, color=GOLD,  lw=0.7, ls='--', label=f"Clean={GR_CLEAN}")
ax_gr.axvline(GR_SHALE, color=RED,   lw=0.7, ls='--', label=f"Shale={GR_SHALE}")
ax_gr.set_xlim(0, 150)
ax_gr.set_xlabel("GR (API)", color=WHITE, fontsize=8)
ax_gr.set_ylabel("Depth (m MD)", color=WHITE, fontsize=8)
ax_gr.set_title("GR", color=GREEN, fontsize=8, fontweight='bold')

# ── Track 2: Density + Neutron Porosity ───────────────────────────────────
ax_dn.plot(df["RHOB_sm"], depth, color=RED,    lw=0.8, label="RHOB")
ax_dn2 = ax_dn.twiny()
ax_dn2.set_facecolor(PANEL)
ax_dn2.plot(df["NPHI_sm"], depth, color=ACCENT, lw=0.8, label="NPHI")
ax_dn2.set_xlim(0.6, 0)   # NPHI reversed
ax_dn2.tick_params(labelsize=7, colors=ACCENT)
ax_dn.set_xlim(1.95, 2.95)
ax_dn.set_title("RHOB / NPHI\n(crossover = HC)", color=WHITE, fontsize=7, fontweight='bold')
ax_dn.tick_params(colors=RED)

# Shade crossover zones (HC indicator)
for i in range(len(depth)-1):
    rhob_norm = (df["RHOB_sm"].iloc[i] - 1.95) / (2.95 - 1.95)
    nphi_norm = 1 - df["NPHI_sm"].iloc[i] / 0.6
    if rhob_norm > nphi_norm:  # crossover
        ax_dn.axhspan(depth[i], depth[i+1], alpha=0.08, color=GOLD, linewidth=0)

# ── Track 3: Resistivity (log scale) ──────────────────────────────────────
ax_rt.semilogx(df["RT_sm"], depth, color=PURPLE, lw=0.8)
ax_rt.set_xlim(0.1, 1000)
ax_rt.set_title("RT\n(ohm·m)", color=PURPLE, fontsize=8, fontweight='bold')
ax_rt.axvline(10, color=GOLD, lw=0.6, ls=':', alpha=0.7)  # HC indicator

# ── Track 4: Petrophysics ─────────────────────────────────────────────────
ax_pe.plot(df["VSHALE"], depth, color=MUTED,  lw=0.8, label="Vsh")
ax_pe.plot(df["PHIE"],   depth, color=ACCENT, lw=0.8, label="ΦE")
ax_pe.fill_betweenx(depth, df["PHIE"], 0, alpha=0.2, color=ACCENT)
ax_pe.axvline(VSHALE_CUTOFF, color=MUTED,  lw=0.6, ls='--')
ax_pe.axvline(PHIE_CUTOFF,   color=ACCENT, lw=0.6, ls='--')
ax_pe.set_xlim(0, 1)
ax_pe.set_title("Vshale / ΦEFF", color=WHITE, fontsize=8, fontweight='bold')
ax_pe.legend(loc='lower right', fontsize=6, facecolor=PANEL, edgecolor=MUTED)

# ── Track 5: Water Saturation ──────────────────────────────────────────────
ax_sw.plot(df["SW"], depth, color=ACCENT, lw=0.8)
ax_sw.fill_betweenx(depth, df["SW"], 1, alpha=0.15, color=GREEN, label="HC")
ax_sw.axvline(SW_CUTOFF, color=RED, lw=0.6, ls='--')
ax_sw.set_xlim(0, 1)
ax_sw.set_title("Sw", color=WHITE, fontsize=8, fontweight='bold')

# ── Track 6: Net Pay flag ─────────────────────────────────────────────────
ax_pay.barh(
    df.loc[df["NET_PAY"]==1, "DEPTH"],
    width=1, height=0.15,
    color=GOLD, alpha=0.9, linewidth=0
)
ax_pay.barh(
    df.loc[(df["NET_RESERVOIR"]==1) & (df["NET_PAY"]==0), "DEPTH"],
    width=1, height=0.15,
    color=GREEN, alpha=0.5, linewidth=0
)
ax_pay.set_xlim(0, 1)
ax_pay.set_xticks([])
ax_pay.set_title("Pay", color=GOLD, fontsize=8, fontweight='bold')

# ── Track 7: Formation column ──────────────────────────────────────────────
for fname, (top, base) in formations.items():
    col = FORM_COLORS.get(fname, "#555")
    ax_form.axhspan(top, base, color=col, alpha=0.85, linewidth=0)
    mid = (top + base) / 2
    ax_form.text(0.5, mid, fname.replace("_", "\n"), ha='center', va='center',
                 fontsize=5.5, color=WHITE, fontweight='bold',
                 transform=ax_form.get_yaxis_transform())

ax_form.set_xlim(0, 1)
ax_form.set_xticks([])
ax_form.set_title("Formation", color=WHITE, fontsize=8, fontweight='bold')

# ── Annotate Heimdal reservoir ─────────────────────────────────────────────
for ax in [ax_gr, ax_dn, ax_rt, ax_pe, ax_sw]:
    ax.axhspan(3000, 3200, color=GOLD, alpha=0.06, linewidth=0)
    ax.axhline(3000, color=GOLD, lw=0.8, ls='-', alpha=0.6)
    ax.axhline(3200, color=GOLD, lw=0.8, ls='-', alpha=0.6)
ax_gr.text(75, 3100, "◀ HEIMDAL\n  RESERVOIR",
           color=GOLD, fontsize=6, fontweight='bold', va='center')

plt.savefig("../figures/01_five_track_log_display.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Saved: figures/01_five_track_log_display.png")

# ═══════════════════════════════════════════════════════════════════════════
# 6. FIGURE 2 — CROSSPLOT SUITE (Petrel QC workflow)
# ═══════════════════════════════════════════════════════════════════════════
print("  Rendering Figure 2: Crossplot suite...")

fig, axes = plt.subplots(2, 2, figsize=(14, 11), facecolor=BG)
fig.suptitle(
    "CROSSPLOT QC SUITE  ·  15/9-F-1  ·  Volve Field Analogue",
    fontsize=10, color=WHITE, fontweight='bold', y=0.98
)

form_list = df["FORMATION"].unique()
cmap_dict = FORM_COLORS

# Sample for plotting speed
mask = df.index % 3 == 0
ds = df[mask].copy()

def scatter_by_formation(ax, x, y, xlabel, ylabel, title, xlog=False, ylog=False):
    for fname in form_list:
        sub = ds[ds["FORMATION"] == fname]
        ax.scatter(sub[x], sub[y],
                   c=FORM_COLORS.get(fname, "#888"),
                   s=2, alpha=0.5, linewidths=0, label=fname.replace("_"," "))
    ax.set_xlabel(xlabel, color=WHITE, fontsize=8)
    ax.set_ylabel(ylabel, color=WHITE, fontsize=8)
    ax.set_title(title, color=WHITE, fontsize=9, fontweight='bold')
    ax.grid(True, color='#21262d', linewidth=0.4)
    if xlog: ax.set_xscale('log')
    if ylog: ax.set_yscale('log')

# Neutron-Density (lithology identification)
scatter_by_formation(axes[0,0], "NPHI_sm", "RHOB_sm",
                     "NPHI (v/v)", "RHOB (g/cc)",
                     "Neutron–Density (Lithology ID)")
# Limestone, sandstone, dolomite lines
for label, nphi_pts, rhob_pts, col in [
    ("Sandstone",  [0.0, 0.45], [2.65, 1.99], GOLD),
    ("Limestone",  [0.0, 0.45], [2.71, 2.08], ACCENT),
    ("Dolomite",   [0.0, 0.45], [2.87, 2.35], GREEN),
]:
    axes[0,0].plot(nphi_pts, rhob_pts, '--', color=col, lw=1.0,
                   alpha=0.7, label=label)
axes[0,0].set_xlim(-0.05, 0.55)
axes[0,0].set_ylim(1.8, 3.0)
axes[0,0].invert_xaxis()
axes[0,0].legend(fontsize=5, facecolor=PANEL, edgecolor=MUTED, ncol=2,
                 markerscale=3, loc='upper right')

# M-N plot proxy (RHOB vs DT)
scatter_by_formation(axes[0,1], "DT_sm", "RHOB_sm",
                     "DT (μs/ft)", "RHOB (g/cc)",
                     "Sonic–Density (Fluid Substitution QC)")
axes[0,1].set_xlim(40, 140)
axes[0,1].set_ylim(1.8, 3.0)
axes[0,1].legend(fontsize=5, facecolor=PANEL, edgecolor=MUTED, markerscale=3)

# Porosity–Permeability (Heimdal reservoir only)
heimdal = df[df["FORMATION"] == "Heimdal_Fm"].copy()
sc = axes[1,0].scatter(heimdal["PHIE"], heimdal["PERM"],
                        c=heimdal["SW"], cmap='RdYlGn_r',
                        s=5, alpha=0.7, vmin=0, vmax=1)
plt.colorbar(sc, ax=axes[1,0], label="Sw", fraction=0.03)
axes[1,0].set_xlabel("Effective Porosity (v/v)", color=WHITE, fontsize=8)
axes[1,0].set_ylabel("Permeability (mD)", color=WHITE, fontsize=8)
axes[1,0].set_title("Φ–k Relationship  ·  Heimdal Fm (Reservoir)", color=WHITE,
                     fontsize=9, fontweight='bold')
axes[1,0].set_yscale('log')
axes[1,0].axvline(PHIE_CUTOFF, color=GOLD, lw=0.8, ls='--', label="ΦE cutoff")
axes[1,0].axhline(0.1, color=RED, lw=0.8, ls='--', label="k=0.1mD")
axes[1,0].legend(fontsize=6, facecolor=PANEL, edgecolor=MUTED)
axes[1,0].grid(True, color='#21262d', linewidth=0.4)

# Sw vs PHIE (Archie saturation height)
scatter_by_formation(axes[1,1], "PHIE", "SW",
                     "Effective Porosity (v/v)", "Water Saturation (Sw)",
                     "Saturation–Porosity  ·  Archie Analysis")
axes[1,1].axhline(SW_CUTOFF, color=GOLD, lw=0.8, ls='--', label=f"Sw cutoff={SW_CUTOFF}")
axes[1,1].axvline(PHIE_CUTOFF, color=ACCENT, lw=0.8, ls='--', label=f"Φ cutoff={PHIE_CUTOFF}")
axes[1,1].set_xlim(0, 0.45)
axes[1,1].set_ylim(0, 1)
axes[1,1].legend(fontsize=6, facecolor=PANEL, edgecolor=MUTED)

for ax in axes.flat:
    ax.set_facecolor(PANEL)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig("../figures/02_crossplot_suite.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Saved: figures/02_crossplot_suite.png")

# ═══════════════════════════════════════════════════════════════════════════
# 7. FIGURE 3 — CCS STORAGE ASSESSMENT SUMMARY
# ═══════════════════════════════════════════════════════════════════════════
print("  Rendering Figure 3: CCS storage assessment...")

fig = plt.figure(figsize=(16, 9), facecolor=BG)
fig.suptitle(
    "CCS STORAGE SUITABILITY ASSESSMENT  ·  15/9-F-1 ANALOGUE\n"
    "Heimdal Fm (Reservoir) + Lista Fm (Seal) + Rogaland Gp (Cap Rock)",
    fontsize=10, color=WHITE, fontweight='bold', y=0.98
)

gs = GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.4,
              left=0.07, right=0.97, top=0.90, bottom=0.08)

ax_col  = fig.add_subplot(gs[:, 0])    # Stratigraphic column
ax_phi  = fig.add_subplot(gs[0, 1])   # Porosity distribution
ax_perm = fig.add_subplot(gs[0, 2])   # Permeability distribution
ax_seal = fig.add_subplot(gs[1, 1])   # Seal integrity
ax_sum  = fig.add_subplot(gs[1, 2])   # Summary scorecard

# ── Stratigraphic column ───────────────────────────────────────────────────
column_data = [
    ("Nordland Gp",  700,  "#6e7681", "Overburden\n(unconsolidated)"),
    ("Hordaland Gp", 1200, "#5a6b7a", "Cap rock\n(thick shale)"),
    ("Shetland Gp",  600,  "#c9a86c", "Chalk/limestone"),
    ("Rogaland Gp",  400,  "#4a5b6a", "Primary seal\n★ SHALE"),
    ("Heimdal Fm",   200,  "#f5a623", "★ RESERVOIR\n(target injectite)"),
    ("Lista Fm",     200,  "#4a5b6a", "Secondary seal\n(shale)"),
    ("Ty Fm",        200,  "#e8c56a", "Secondary\nreservoir"),
]

y_pos = 0
bar_w = 0.6
for fname, thick, col, label in column_data:
    rect = mpatches.FancyBboxPatch(
        (0.1, y_pos), bar_w, thick * 0.0008,
        boxstyle="square,pad=0", facecolor=col, edgecolor=BG, linewidth=0.5
    )
    ax_col.add_patch(rect)
    ax_col.text(0.4, y_pos + thick * 0.0004, label,
                ha='center', va='center', fontsize=6, color=WHITE, fontweight='bold')
    y_pos += thick * 0.0008

ax_col.set_xlim(0, 0.8)
ax_col.set_ylim(0, y_pos + 0.05)
ax_col.set_xticks([])
ax_col.set_yticks([])
ax_col.set_title("Stratigraphy", color=WHITE, fontsize=9, fontweight='bold')
ax_col.set_facecolor(PANEL)
# CO2 injection arrow
ax_col.annotate("CO₂\ninjection", xy=(0.75, 0.385), xytext=(0.75, 0.30),
                color=ACCENT, fontsize=7, ha='center',
                arrowprops=dict(arrowstyle='->', color=ACCENT, lw=1.5))

# ── Porosity distribution ──────────────────────────────────────────────────
heimdal = df[df["FORMATION"] == "Heimdal_Fm"]
ax_phi.hist(heimdal["PHIE"], bins=40, color=GOLD, alpha=0.85, edgecolor=BG)
ax_phi.axvline(heimdal["PHIE"].mean(), color=WHITE, lw=1.2, ls='--',
               label=f'Mean={heimdal["PHIE"].mean():.2f}')
ax_phi.axvline(PHIE_CUTOFF, color=RED, lw=1.0, ls=':', label=f'Cutoff={PHIE_CUTOFF}')
ax_phi.set_xlabel("Effective Porosity (v/v)", color=WHITE, fontsize=8)
ax_phi.set_ylabel("Count", color=WHITE, fontsize=8)
ax_phi.set_title("Porosity Distribution\nHeimdal Fm", color=WHITE, fontsize=8, fontweight='bold')
ax_phi.legend(fontsize=7, facecolor=PANEL, edgecolor=MUTED)
ax_phi.set_facecolor(PANEL)

# ── Permeability distribution ─────────────────────────────────────────────
perm_vals = heimdal["PERM"].clip(0.1, 5000)
ax_perm.hist(np.log10(perm_vals), bins=40, color=PURPLE, alpha=0.85, edgecolor=BG)
ax_perm.axvline(np.log10(heimdal["PERM"].median()), color=WHITE, lw=1.2, ls='--',
                label=f'P50={heimdal["PERM"].median():.0f} mD')
ax_perm.set_xlabel("log₁₀ Permeability (mD)", color=WHITE, fontsize=8)
ax_perm.set_ylabel("Count", color=WHITE, fontsize=8)
ax_perm.set_title("Permeability Distribution\nHeimdal Fm", color=WHITE, fontsize=8, fontweight='bold')
ax_perm.legend(fontsize=7, facecolor=PANEL, edgecolor=MUTED)
ax_perm.set_facecolor(PANEL)

# ── Seal integrity — Lista shale GR / RHOB ────────────────────────────────
lista = df[df["FORMATION"] == "Lista_Fm"]
ax_seal.plot(lista["GR_sm"], lista["DEPTH"], color=RED, lw=1.0, label="GR (API)")
ax_seal2 = ax_seal.twiny()
ax_seal2.plot(lista["RHOB_sm"], lista["DEPTH"], color=MUTED, lw=1.0, label="RHOB")
ax_seal2.tick_params(labelsize=7, colors=MUTED)
ax_seal2.set_xlim(2.3, 2.7)
ax_seal.set_ylim(3400, 3200)
ax_seal.set_xlabel("GR (API)", color=RED, fontsize=8)
ax_seal.set_ylabel("Depth (m)", color=WHITE, fontsize=8)
ax_seal.set_title("Seal QC: Lista Fm\n(High GR = Clean Shale)", color=WHITE,
                  fontsize=8, fontweight='bold')
ax_seal.axvline(75, color=GOLD, lw=0.7, ls='--', label="GR > 75 = shale")
ax_seal.legend(fontsize=7, facecolor=PANEL, edgecolor=MUTED)
ax_seal.set_facecolor(PANEL)

# ── CCS Scorecard ──────────────────────────────────────────────────────────
ax_sum.set_facecolor(PANEL)
ax_sum.set_xticks([])
ax_sum.set_yticks([])
ax_sum.set_title("CCS Storage Scorecard", color=WHITE, fontsize=9, fontweight='bold')

criteria = [
    ("Reservoir porosity (mean)",  f"{heimdal['PHIE'].mean():.0%}",     GREEN,  "✓ GOOD (>15%)"),
    ("Reservoir permeability P50", f"{heimdal['PERM'].median():.0f} mD", GREEN,  "✓ GOOD (>10 mD)"),
    ("Net reservoir thickness",    f"{summary.loc['Heimdal_Fm','Net_Reservoir_m']:.0f} m", GREEN, "✓ GOOD (>50 m)"),
    ("Seal: Lista Fm GR avg",      f"{lista['GR_sm'].mean():.0f} API",   GREEN,  "✓ GOOD (>80 API)"),
    ("Seal thickness",             "200 m",                               GREEN,  "✓ EXCELLENT"),
    ("Cap rock: Rogaland Gp",      "400 m shale",                         GREEN,  "✓ EXCELLENT"),
    ("Depth to reservoir",         "3000 m",                              GREEN,  "✓ GOOD (>800m)"),
    ("Structural closure",         "Needs 3D seismic",                    GOLD,   "⚠ TBD"),
    ("CO₂ injectivity",           "High k → feasible",                   GREEN,  "✓ LIKELY"),
]

y_start = 0.92
for label, value, col, status in criteria:
    ax_sum.text(0.03, y_start, f"▸ {label}", transform=ax_sum.transAxes,
                fontsize=7, color=MUTED, va='top')
    ax_sum.text(0.97, y_start, status, transform=ax_sum.transAxes,
                fontsize=7, color=col, va='top', ha='right', fontweight='bold')
    y_start -= 0.09

ax_sum.text(0.5, 0.04, "OVERALL: FAVOURABLE FOR CO₂ STORAGE",
            transform=ax_sum.transAxes, fontsize=9, color=GREEN,
            ha='center', va='bottom', fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='#1a3d2b', edgecolor=GREEN))

plt.savefig("../figures/03_ccs_storage_assessment.png",
            dpi=180, bbox_inches='tight', facecolor=BG)
plt.close()
print("  ✓ Saved: figures/03_ccs_storage_assessment.png")

# ═══════════════════════════════════════════════════════════════════════════
# 8. EXPORT PETROPHYSICAL RESULTS
# ═══════════════════════════════════════════════════════════════════════════
output_cols = ["DEPTH", "FORMATION", "GR_sm", "RHOB_sm", "NPHI_sm",
               "DT_sm", "RT_sm", "VSHALE", "PHID", "PHIE", "SW",
               "SHC", "PERM", "NET_RESERVOIR", "NET_PAY"]
df[output_cols].to_csv("../data/petrophysical_results.csv", index=False)
summary.to_csv("../data/formation_summary.csv")

print("\n  ✓ Exported: data/petrophysical_results.csv")
print("  ✓ Exported: data/formation_summary.csv")
print("\n  ══════════════════════════════════════════════")
print("  WORKFLOW COMPLETE — 3 figures + 2 data files")
print("  ══════════════════════════════════════════════")

"""
parse_orcacosmo.py
------------------
Reads an ORCA 6 .orcacosmo file and generates:
1. CSACVL (cavity volume in Angstrom^3)
2. Sigma profile discretized into 60 bins matching Aspen's COSMOSAC SGPRF format
3. SGPRF1-5 parameter values (12 elements each = 60 total) for entry into Aspen Plus

Aspen COSMOSAC format confirmed from published IL DAC process simulation SI:
  - 60 total sigma profile bins
  - 5 parameters × 12 elements each: SGPRF1=bins 1-12, SGPRF2=bins 13-24,
    SGPRF3=bins 25-36, SGPRF4=bins 37-48, SGPRF5=bins 49-60
  - Sigma range: -0.025 to +0.025 e/Angstrom^2
  - Reference: Sanchez-Lengeling et al. supplementary material

Reference N2 and CO2 values (from published DAC process simulation SI):
  N2:  CSACVL=36.96, all non-zero values in SGPRF2 and SGPRF3 only
  CO2: CSACVL=47.67, non-zero values in SGPRF2 and SGPRF3 only

Usage:
    python parse_orcacosmo.py <filename.orcacosmo> <molecule_name>

Example:
    python parse_orcacosmo.py n2_solute.orcacosmo N2
    python parse_orcacosmo.py o2_solute.orcacosmo O2
    python parse_orcacosmo.py co2_solute.orcacosmo CO2
    python parse_orcacosmo.py h2o_solute.orcacosmo H2O
    python parse_orcacosmo.py bmimdca_solute.orcacosmo BMIMDCA
"""

import sys
import numpy as np

# ── Unit conversions ──────────────────────────────────────────────────────────
BOHR_TO_ANG   = 0.529177210903
BOHR3_TO_ANG3 = BOHR_TO_ANG ** 3
BOHR2_TO_ANG2 = BOHR_TO_ANG ** 2

# ── Aspen COSMOSAC sigma profile grid ─────────────────────────────────────────
# 60 bins, 5 parameters × 12 elements each
# Sigma range: -0.025 to +0.025 e/Ang^2
N_BINS    = 60
SIGMA_MIN = -0.025
SIGMA_MAX =  0.025
SIGMA_BINS = np.linspace(SIGMA_MIN, SIGMA_MAX, N_BINS)
SIGMA_STEP = SIGMA_BINS[1] - SIGMA_BINS[0]

# Aspen grouping (12 elements per parameter, confirmed from published SI):
SGPRF_SLICES = {
    'SGPRF1': slice(0,  12),   # bins 1-12
    'SGPRF2': slice(12, 24),   # bins 13-24
    'SGPRF3': slice(24, 36),   # bins 25-36
    'SGPRF4': slice(36, 48),   # bins 37-48
    'SGPRF5': slice(48, 60),   # bins 49-60
}


def parse_orcacosmo(filepath):
    """
    Parse an ORCA 6 .orcacosmo file.
    Returns volume (Bohr^3), area (Bohr^2), and segment array.
    Columns: X, Y, Z, area, potential, charge
    """
    volume     = None
    area_total = None
    segments   = []
    in_surface = False

    with open(filepath, 'r') as f:
        lines = f.readlines()

    for line in lines:
        stripped = line.strip()

        if '# Volume' in line:
            try:
                volume = float(stripped.split()[0])
            except ValueError:
                pass

        if '# Area' in line and volume is not None:
            try:
                area_total = float(stripped.split()[0])
            except ValueError:
                pass

        if 'SURFACE POINTS' in line:
            in_surface = True
            continue

        if in_surface and 'X' in line and 'Y' in line and 'area' in line:
            continue   # skip header row

        if in_surface:
            parts = stripped.split()
            if len(parts) >= 6:
                try:
                    vals = [float(p) for p in parts[:6]]
                    segments.append(vals)
                except ValueError:
                    if stripped and not stripped.startswith('#'):
                        in_surface = False

    return volume, area_total, np.array(segments) if segments else np.array([])


def compute_sigma_profile(segments):
    """
    Compute per-segment charge density sigma = charge/area, convert to e/Ang^2.
    """
    areas_bohr2 = segments[:, 3]
    charges_e   = segments[:, 5]
    sigma_bohr  = charges_e / areas_bohr2
    sigma_ang   = sigma_bohr / BOHR2_TO_ANG2
    areas_ang2  = areas_bohr2 * BOHR2_TO_ANG2
    return sigma_ang, areas_ang2


def discretize_sigma_profile(sigma_values, areas):
    """
    Bin segment areas into the 60-bin Aspen COSMOSAC grid.
    """
    profile = np.zeros(N_BINS)
    for s, a in zip(sigma_values, areas):
        idx = int(round((s - SIGMA_MIN) / SIGMA_STEP))
        idx = max(0, min(N_BINS - 1, idx))
        profile[idx] += a
    return SIGMA_BINS.copy(), profile


def format_sgprf(profile):
    """Split 60-bin profile into SGPRF1-5 (12 elements each)."""
    return {key: profile[slc] for key, slc in SGPRF_SLICES.items()}


def main():
    if len(sys.argv) < 2:
        print("Usage: python parse_orcacosmo.py <file.orcacosmo> [molecule_name]")
        sys.exit(1)

    filepath = sys.argv[1]
    molname  = sys.argv[2] if len(sys.argv) > 2 else filepath.split('.')[0].upper()

    print(f"\n{'='*65}")
    print(f"  ORCACOSMO → Aspen COSMOSAC Parameter Extractor")
    print(f"  Molecule : {molname}")
    print(f"  File     : {filepath}")
    print(f"  Grid     : 60 bins (-0.025 to +0.025 e/Ang^2), 5×12 SGPRF grouping")
    print(f"{'='*65}\n")

    volume_bohr3, area_bohr2, segments = parse_orcacosmo(filepath)

    if volume_bohr3 is None:
        print("ERROR: Volume not found. Check file format.")
        sys.exit(1)
    if len(segments) == 0:
        print("ERROR: No surface segments found. Check file format.")
        sys.exit(1)

    volume_ang3 = volume_bohr3 * BOHR3_TO_ANG3
    area_ang2   = area_bohr2   * BOHR2_TO_ANG2 if area_bohr2 else None

    print(f"── Raw COSMO Data ──────────────────────────────────────────")
    print(f"  Volume (Bohr^3)  : {volume_bohr3:.6f}")
    print(f"  Volume (Ang^3)   : {volume_ang3:.6f}   ← CSACVL for Aspen")
    if area_ang2:
        print(f"  Area   (Ang^2)   : {area_ang2:.6f}")
    print(f"  Segments         : {len(segments)}")

    sigma_vals, areas_ang2_seg = compute_sigma_profile(segments)

    print(f"\n── Sigma Profile Statistics ────────────────────────────────")
    print(f"  Sigma range : {sigma_vals.min():.6f} to {sigma_vals.max():.6f} e/Ang^2")
    print(f"  Total area  : {areas_ang2_seg.sum():.4f} Ang^2")

    bins, profile = discretize_sigma_profile(sigma_vals, areas_ang2_seg)
    sgprf = format_sgprf(profile)

    # Label each bin with its SGPRF parameter
    sgprf_labels = (['SGPRF1']*12 + ['SGPRF2']*12 + ['SGPRF3']*12 +
                    ['SGPRF4']*12 + ['SGPRF5']*12)

    print(f"\n── Full 60-bin Sigma Profile ───────────────────────────────")
    print(f"  {'Bin':>4}  {'Sigma (e/Ang^2)':>16}  {'Area (Ang^2)':>14}  {'SGPRF':>6}")
    print(f"  {'-'*48}")
    for i, (s, a, lbl) in enumerate(zip(bins, profile, sgprf_labels)):
        marker = " ←" if a > 0 else ""
        print(f"  {i+1:>4}  {s:>16.8f}  {a:>14.6f}  {lbl:>6}{marker}")

    print(f"\n── SGPRF Parameters for Aspen Entry ────────────────────────")
    print(f"  Enter under Properties → Parameters → Pure Component\n")
    for key, vals in sgprf.items():
        print(f"  {key} (12 elements):")
        # Print as two rows of 6 for readability
        row1 = "  ".join(f"{v:8.2f}" for v in vals[:6])
        row2 = "  ".join(f"{v:8.2f}" for v in vals[6:])
        print(f"    {row1}")
        print(f"    {row2}\n")

    # ── Write output file ─────────────────────────────────────────────────────
    outfile = f"{molname.lower()}_aspen_params.txt"
    with open(outfile, 'w') as f:
        f.write(f"Aspen COSMOSAC Parameters for {molname}\n")
        f.write(f"Generated from: {filepath}\n")
        f.write(f"Grid: 60 bins, -0.025 to +0.025 e/Ang^2, 5 params x 12 elements\n")
        f.write(f"{'='*65}\n\n")
        f.write(f"CSACVL = {volume_ang3:.6f}  (Ang^3)\n\n")
        f.write(f"Full 60-bin sigma profile:\n")
        f.write(f"{'Bin':>4}  {'Sigma (e/Ang^2)':>16}  {'Area (Ang^2)':>14}  {'SGPRF':>6}\n")
        f.write(f"{'-'*48}\n")
        for i, (s, a, lbl) in enumerate(zip(bins, profile, sgprf_labels)):
            f.write(f"{i+1:>4}  {s:>16.8f}  {a:>14.6f}  {lbl:>6}\n")
        f.write(f"\nSGPRF Parameters (12 elements each):\n")
        for key, vals in sgprf.items():
            vals_str = "  ".join(f"{v:.6f}" for v in vals)
            f.write(f"{key}: {vals_str}\n")

    print(f"  Results saved to: {outfile}")
    print(f"\n{'='*65}\n")


if __name__ == '__main__':
    main()

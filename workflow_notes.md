# Workflow Notes — Petrophysical Analysis Decisions

## Log QC Strategy
- 5-point uniform smoothing applied to all logs before petrophysical computation
- Removes high-frequency noise from borehole rugosity and tool noise
- Equivalent to Petrel's "Well Log Filter" function

## Vshale Method Selection
Larionov (1969) non-linear equation selected for Tertiary-age shales:
- Vsh = 0.083 × (2^(3.7 × IGR) − 1)
- More appropriate than linear IGR for younger, less compacted shales
- Reduces Vshale overestimation in carbonate-rich intervals

## Porosity Method
- PHID (density porosity) computed as primary porosity indicator
- PHIE (effective porosity) = average of PHID and NPHI — reduces fluid effects
- Formation water density assumed 1.05 g/cc (NaCl brine, ~80°C, ~30,000 ppm)

## Archie Parameters
Standard parameters used for clean consolidated sandstone:
- a=1.0 (tortuosity): appropriate for intergranular porosity
- m=2.0 (cementation): appropriate for well-cemented sandstone
- n=2.0 (saturation): standard assumption, would be refined with special core analysis

## Net Pay Cutoffs — Basis
Cutoffs are conservative estimates appropriate for a CCS storage assessment context:
- Vshale < 0.30: Standard reservoir quality cutoff for clastic reservoirs
- PHIE > 0.10: Minimum porosity for meaningful CO₂ storage and injectivity
- Sw < 0.65: Ensures residual hydrocarbon/gas saturation for storage capacity

## CCS Context
For CCS storage (vs. hydrocarbon production), Sw cutoff is less critical —
CO₂ can be injected regardless of initial fluid saturation. The primary
constraints are porosity, permeability, and seal integrity. The Sw cutoff
here is used to identify original hydrocarbon columns as a proxy for
structural trap integrity and seal quality.

## Known Limitations
1. Single well — no lateral heterogeneity captured
2. Synthetic data — calibrated to Volve stratigraphy but not actual Volve logs
3. No special core analysis (SCAL) data — Archie parameters are assumed
4. No temperature/pressure corrections applied to log readings
5. No geomechanical analysis of seal integrity under CO₂ injection conditions

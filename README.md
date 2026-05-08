# Reservoir Geology Portfolio — Volve Field Analogue
### Well Log Analysis · Petrophysical Characterisation · CCS Storage Assessment

[GitHub](https://github.com/Joexy1286)

---

## Overview

This repository demonstrates a complete **reservoir geologist workflow** applied to a synthetic well log dataset modelled on the publicly available **Equinor Volve field** (Norwegian North Sea). The workflow mirrors industry practice as performed in Petrel, implemented here in Python for reproducibility and open-source accessibility.

The project is framed as a **CCS (CO₂ storage) suitability assessment** — directly reflecting the energy transition context in which reservoir geology skills are increasingly deployed.

---

## Workflow Summary

```
Raw LAS-format logs
       │
       ▼
  1. Log QC & smoothing
       │
       ▼
  2. Lithology identification (GR, neutron-density crossplot)
       │
       ▼
  3. Petrophysical computation
      ├── Vshale (Larionov non-linear)
      ├── Density porosity (PHID)
      ├── Effective porosity (PHIE)
      ├── Water saturation (Archie)
      └── Permeability (Timur-Coates)
       │
       ▼
  4. Net reservoir & net pay flagging
       │
       ▼
  5. Crossplot QC suite
      ├── Neutron-Density (lithology)
      ├── Sonic-Density (fluid substitution)
      ├── Φ–k relationship (reservoir quality)
      └── Sw–Φ (Archie saturation height)
       │
       ▼
  6. CCS storage suitability assessment
      ├── Reservoir: Heimdal Fm (sand injectite, 3000–3200 m)
      ├── Seal: Lista Fm (shale, 3200–3400 m)
      └── Cap rock: Rogaland Gp (400 m thick shale)
```

---

## Key Results

| Formation | Depth (m) | Mean ΦEFF | k P50 (mD) | Net Res. (m) | Role |
|---|---|---|---|---|---|
| **Heimdal Fm** | 3000–3200 | **25%** | **140 mD** | **200 m** | **Primary reservoir/injectite** |
| Lista Fm | 3200–3400 | 15% | 13 mD | — |  Primary seal |
| Rogaland Gp | 2600–3000 | 17% | 27 mD | — |  Cap rock |
| Ty Fm | 3400–3600 | 22% | 76 mD | 200 m | Secondary reservoir |

### CCS Assessment: Heimdal Fm → **FAVOURABLE**
-  Porosity 25% (target >15%)
-  Permeability P50 = 140 mD (target >10 mD for injectivity)
-  Net reservoir 200 m (target >50 m)
-  Seal thickness 200 m
-  Depth 3000 m (supercritical CO₂ conditions assured above ~800 m)
-  Structural closure requires 3D seismic confirmation

---

## Figures

### Figure 1 — Five-Track Log Display
*GR · Density-Neutron · Resistivity · Petrophysics · Net Pay*

![Five-track log display](figures/01_five_track_log_display.png)

**Reading the log:**
- **Track 1 (GR, green):** Low GR in Heimdal Fm (3000–3200 m) confirms clean sand
- **Track 2 (Density-Neutron):** Crossover between RHOB and NPHI in pay zone indicates hydrocarbon effect
- **Track 3 (Resistivity, log scale):** Elevated RT >10 ohm·m in Heimdal confirms hydrocarbon saturation
- **Track 4 (Vshale/ΦEFF):** Effective porosity peaks at ~25% in reservoir interval
- **Track 6 (Net Pay, gold):** Continuous net pay through Heimdal Fm

---

### Figure 2 — Crossplot QC Suite

![Crossplot suite](figures/02_crossplot_suite.png)

Four standard Petrel QC crossplots:
- **Neutron-Density (top left):** Points colour-coded by formation. Heimdal Fm plots near the sandstone line, confirming clean quartz sand. Crossover cluster above the sandstone line indicates gas effect.
- **Sonic-Density (top right):** Used for fluid substitution and rock physics QC. Heimdal sand cluster shows characteristic low-density/high-DT signature.
- **Φ–k (bottom left):** Heimdal Fm only, colour-coded by Sw. Strong porosity-permeability correlation with good injectivity in the low-Sw (hydrocarbon-bearing) zone.
- **Saturation–Porosity (bottom right):** Archie relationship. Heimdal sand clearly separates from shales with low Sw and high porosity — confirming net pay classification.

---

### Figure 3 — CCS Storage Assessment

![CCS storage assessment](figures/03_ccs_storage_assessment.png)

Integrated storage suitability assessment using:
- Stratigraphic column showing reservoir-seal-caprock trapping geometry
- Porosity and permeability distributions for injectivity modelling
- Seal QC using Lista Fm shale indicators
- Summary scorecard against IEA/IEAGHG storage criteria

---

## Methods & Parameters

### Petrophysical Parameters
| Parameter | Value | Basis |
|---|---|---|
| GR clean baseline | 20 API | Heimdal Fm clean sand |
| GR shale baseline | 110 API | Rogaland/Lista shales |
| Vshale equation | Larionov (Tertiary) | Appropriate for Cenozoic shales |
| Matrix density | 2.65 g/cc | Quartz sandstone |
| Fluid density | 1.05 g/cc | Saline formation water |
| Archie: a, m, n | 1.0, 2.0, 2.0 | Standard for consolidated sand |
| Formation water resistivity | 0.04 ohm·m | ~80°C formation temperature |
| Permeability model | Timur-Coates | Appropriate for clastics |

### Net Pay Cutoffs
| Parameter | Cutoff | Basis |
|---|---|---|
| Vshale | < 0.30 | Reservoir quality cutoff |
| ΦEFF | > 0.10 | Minimum storage porosity |
| Sw | < 0.65 | Hydrocarbon saturation cutoff |

---

## Repository Structure

```
reservoir-geology-portfolio/
├── notebooks/
│   ├── 01_well_log_analysis.py      # Main workflow script
│   └── 02_seismic_interpretation.py # (in development)
├── data/
│   ├── generate_synthetic_logs.py   # Synthetic Volve-analogue data generator
│   ├── petrophysical_results.csv    # Computed petrophysical outputs
│   └── formation_summary.csv        # Formation-level statistics
├── figures/
│   ├── 01_five_track_log_display.png
│   ├── 02_crossplot_suite.png
│   └── 03_ccs_storage_assessment.png
└── docs/
    └── workflow_notes.md
```

---

## Dependencies

```bash
pip install numpy pandas matplotlib scipy
```

No proprietary software required. All outputs reproducible from a standard scientific Python environment.

---

## Petrel Course Context

The workflow implemented here follows the sequence taught in the **Seismic Interpretation Using Petrel** course. Key Petrel workflows replicated in Python:
- Log smoothing (equivalent to Petrel's well log filter)
- Neutron-Density crossplot (Petrel: Well Correlation window)
- Phi-k crossplot with Sw colour coding (Petrel: Crossplot tool)
- Net reservoir / net pay computation (Petrel: Petrophysics plugin)
- Formation evaluation and log track display (Petrel: Well Section window)

---

## Connection to Thesis Research

This reservoir geology work directly complements my MSc thesis on **CO₂ mineralisation in basalt systems** (Fraunhofer IEG, PHREEQC/Reaktoro reactive transport modelling). Both workflows address different scales of the same problem:

- **Thesis:** Geochemical fate of injected CO₂ at pore scale (mineral trapping, reaction kinetics)
- **This portfolio:** Reservoir-scale characterisation for storage site selection (structural trapping, capacity, injectivity)

Together they represent the full CCS subsurface workflow: site selection → injection engineering → long-term geochemical monitoring.

---

## Related Repositories

- [co2-mineralisation-phreeqc](https://github.com/Joexy1286/co2-mineralisation-phreeqc) — PHREEQC reactive transport models for CO₂ mineralisation in Icelandic basalt, Diabase Harz, and Harzburgite (thesis work)

---

*Data in this repository is synthetic, generated to replicate the stratigraphy and log character of the publicly available Equinor Volve dataset. No proprietary data is included.*

# Synaptopodin enables directional mechanoadaptation of integrin-based adhesions — simulation code

Computational model of cell detachment under perpendicular (normal) load for
Qu *et al.*, *Synaptopodin enables directional mechanoadaptation of integrin-based adhesions*.

The model represents a cell as an elliptical array of actin stress fibers that
load integrin bonds. Each bond fails through a **catch–slip off-rate** in which
in-plane shear holds bonds in a protected regime while concentrated normal load
accelerates rupture. A single synaptopodin-dependent parameter **C** sets the
sensitivity of the off-rate to concentrated normal load (WT `C = 1.5`,
synaptopodin-deficient `C = 2.7`). From this the model reproduces the
genotype-dependent detachment, the peripheral redistribution of bonds, and the
load-induced change in cell shape reported in the paper.

A full description of the equations is in the manuscript **Methods**
("Computational model of cell detachment", Equations 7–16) and **Supplementary
Note 1** (sensitivity and robustness analyses). All model parameters, values,
and sources are listed in **Supplementary Table 1**.

---

## 1. System requirements

**Software dependencies**

| Component | Used for | Version tested |
|---|---|---|
| MATLAB | Figure 5b–d (`.m` scripts) | R2020b or later (no toolboxes required for the canonical scripts) |
| Python | Sensitivity analysis (`.py` scripts) | 3.9 – 3.12 (tested on 3.12.2) |
| numpy | numerics | 2.4.4 |
| pandas | tables / CSV I/O | 3.0.2 |
| matplotlib | plotting | 3.10.0 |
| openpyxl | reading/writing `.xlsx` (optional) | any recent |

**Operating systems:** platform-independent (pure MATLAB / pure Python).
Tested on macOS (Darwin 25.x). Expected to run unchanged on Linux and Windows.

**Non-standard hardware:** none. Runs on a normal desktop/laptop; no GPU and no
cluster required. The canonical MATLAB figure scripts do **not** require the
Parallel Computing Toolbox (the exploratory `parameter_sweep*.m` drivers do).

## 2. Installation guide

```bash
git clone https://github.com/kenny003003/Synaptopodin-enables-directional-mechanoadaptation-of-integrin-based-adhesions.git
cd Synaptopodin-enables-directional-mechanoadaptation-of-integrin-based-adhesions
pip install numpy pandas matplotlib openpyxl      # Python dependencies
```

MATLAB requires no installation beyond MATLAB itself. **Typical install time**
on a normal desktop: Python dependencies ≈ 1–2 min; MATLAB already installed.

## 3. Demo

A fast demo that needs no long simulation: the sensitivity results are provided
as `results_local_OAT.csv` and `results_sobol.csv`, and the figure is rendered
from them.

```bash
python3 make_sensitivity_figure.py
```

**Expected output:** `sensitivity_figure.png` — a two-panel bar chart. Panel (a)
local one-at-a-time sensitivity, panel (b) Sobol indices. The baseline off-rate
`k_off0` and the shear-protection Hill exponent `n` are the two most influential
parameters; the pre-tension `F0`, linkage stiffness `k0`, and slip-amplification
scale `K0` are negligible. **Expected run time:** a few seconds.

To regenerate the underlying numbers from scratch:

```bash
python3 sensitivity_analysis.py        # writes results_local_OAT.csv, results_sobol.csv
```

**Expected output** (operating point where the knockout reaches ≈50% detachment):
`k_off0` ≈ +9.7 and `n` ≈ −8.0 as the largest local elasticities, and Sobol
total-effect indices `k_off0` ≈ 0.60, `n` ≈ 0.47 with all others ≤ 0.08.
**Expected run time:** ≈ 5–8 min on a normal desktop.

## 4. Instructions for use

This is a Monte-Carlo simulation; it generates its own stress-fiber ensembles
and requires no external input data. To use it, edit the parameter block at the
top of a script (cell geometry, `C`, force levels, replicate counts) and run it.

**Reproduction of the published computational results**

| Script | Reproduces | Approx. run time |
|---|---|---|
| `aspect_ratio_sweep.m` | Fig. 5b — aspect-ratio change under perpendicular load (WT vs KO) | a few minutes |
| `barplot.m` | Fig. 5c — attachment rate vs applied force (WT vs KO) | a few minutes |
| `rose_diagram.m` | Fig. 5d — angular distribution of surviving bonds | a few minutes |
| `sensitivity_analysis.py` | Supplementary Tables 2 & 3 — local + Sobol sensitivity | ≈ 5–8 min |
| `make_sensitivity_figure.py` | Supplementary Fig. 15 | seconds |

Each MATLAB script writes a `.png` of its figure; `aspect_ratio_sweep.m` also
exports the two-curve data to `figure5b_curve_data.xlsx` (a copy is included as
example output). Run times scale with the replicate counts (`N_TRIALS`,
`N_SETS`, `N_CELLS_PER_SET`, `R`) set at the top of each script; reduce them for
a faster check.

## 5. Repository contents

**Canonical scripts (physical parameterization, reproduce the published figures)**
- `aspect_ratio_sweep.m`, `barplot.m`, `rose_diagram.m` — Figure 5b, 5c, 5d
- `sensitivity_analysis.py` — Supplementary Tables 2 & 3
- `make_sensitivity_figure.py` — Supplementary Fig. 15

**Data / example output**
- `results_local_OAT.csv`, `results_sobol.csv` — sensitivity results
- `figure5b_curve_data.xlsx` — Figure 5b curve data

**Exploratory / earlier-iteration scripts** (kept for transparency; *not* used to
produce the published figures, and may use earlier parameter values)
- `parameter_sweep.m`, `parameter_sweep_v2.m`, `parameter_sweep_v3.m` — parameter-exploration drivers (use the Parallel Computing Toolbox)
- `fiber_polarization.m` — fiber-evolution / shape-change visualization
- `equation_2.m` — simplified two-cable analytical force balance (illustrative)
- `simulation.py` — standalone Python implementation of the core model

**Utilities**
- `RosePlot.m` — interactive rose plot from a user-selected Excel file
- `plot_aspect_ratio.m` — plotting helper

## 6. License

Released under the MIT License — see [LICENSE](LICENSE).

## 7. Citation

If you use this code, please cite Qu *et al.*, *Synaptopodin enables directional
mechanoadaptation of integrin-based adhesions* (and this repository).

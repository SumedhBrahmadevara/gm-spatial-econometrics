
# Housing Affordability in the Face of Gentrification — Code

This repository contains the analysis code and minimal data needed to reproduce the core computational steps for the dissertation *Housing Affordability in the Face of Gentrification: An Economic Inquiry* (Greater Manchester, 1991–2021).

> **Reproduction scope**: The repo includes the exact pipeline structure referenced in Appendix C (Stata, Python, R). Large/raw datasets are excluded; small samples are provided for validation.

## Repository structure

```
gm-gentrification-analysis/
├── code/
│   ├── msoa_processing_scripts.do        # Stata: harmonization, merges, SDM estimation
│   ├── gm_tools_improved.py              # Python: spatial controls (greenspace, schools, stations)
│   ├── census_harmonization.R            # R: 1991 ED → 2011 MSOA areal interpolation
│   └── gentrification_analysis.do        # Stata: index construction and model routines
├── data/
│   ├── sample_datasets/                  # small mock samples for tests/CI
│   └── validation_outputs/               # intermediate outputs (small)
├── docs/
│   └── technical_notes.md                # extra implementation notes
├── tests/
│   └── test_python_smoke.py              # minimal CI test for python utilities
├── .github/workflows/ci.yml              # CI: installs python deps and runs tests
├── requirements.txt                      # Python deps
├── R-packages.txt                        # R deps (use renv if preferred)
├── Stata-packages.txt                    # Stata packages used
├── .gitignore
├── LICENSE
└── CITATION.cff
```

## Quickstart

### 1) Python (spatial controls)
```bash
# Create and activate a venv (recommended)
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate

pip install -r requirements.txt

python - <<'PY'
from code.gm_tools_improved import compute_spatial_controls
print("Module imported. (Run with real data in your workflow.)")
PY
```

### 2) R (areal interpolation)
- Install packages listed in `R-packages.txt` (or initialise `renv` and run `renv::restore()`).
- Run: `Rscript code/census_harmonization.R`

### 3) Stata (harmonization / SDM)
- Install packages listed in `Stata-packages.txt`
- Run in Stata:
```
do code/msoa_processing_scripts.do
do code/gentrification_analysis.do
```

## Reproducibility and versions

- Python: see `requirements.txt`
- R: see `R-packages.txt` (or lock via renv)
- Stata: 16+ with `spxtregress`, `sppack`

## License and citation

- License: MIT (see `LICENSE`).
- Cite via `CITATION.cff` or the dissertation.

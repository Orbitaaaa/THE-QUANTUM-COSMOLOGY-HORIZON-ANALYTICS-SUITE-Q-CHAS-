# Reproducibility package (G11-28)

## Environment

* Python 3.13 with numpy, scipy, numba, sympy, matplotlib (versions in `DATA/round11_manifest.json`). `CODE/requirements.txt` lists the
  dependencies.
* Pythia: `pip install pythia8mc==8.317.2` (it ships the PDF grids; their SHA-256 hashes are in the manifest).
* LaTeX (pdflatex) for the papers.

## One command per result

| Result | Command (from `genesis/CODE/scripts`) | Output | Time on 4 cores |
|---|---|---|---|
| Paper A (all physics numbers) | `python round11_physics.py` | `DATA/round11_physics.json`, `FIGURES/round11_physics.png` | ~1.5 h |
| Causal schematic | `python round11_causal.py` | `FIGURES/round11_causal.png` | seconds |
| Greybody tables | `python round11_greybody.py --workers 4` | `DATA/round11_greybody.json` | ~3 min |
| Background and signal samples | `python round10a.py`, then `python round11a.py --workers 3` | `DATA/round10_events/*.npz` (not committed) | several hours |
| Collider budget (Paper B) | `python round11b.py --workers 3` | `DATA/round11_collider.json`, `DATA/round11_manifest.json`, `ROUND11/collider/analysis_config.json`, `FIGURES/round11_collider.png` | ~1 h |
| Timing demonstrator | `python round11_timing.py` | `DATA/round11_timing.json`, `FIGURES/round11_timing.png` | minutes |
| Tests | `cd ..; python -m pytest tests -q` | — | ~minutes (the long tests are marked) |
| Papers | `cd ../../PAPER; pdflatex paperA_inner_horizon.tex; pdflatex paperB_collider.tex` | PDFs | seconds |

## Verifying regenerated samples

`DATA/round11_manifest.json` lists the SHA-256 hash, event count and cross-section of every sample. Bit-identical regeneration needs
the same Pythia version and platform. Otherwise compare cross-sections and the reach numbers, whose tolerance is set by the MC-statistics
nuisance.

## Archiving

* Tag the commit used for each paper.
* Deposit the code, JSON data and event samples in a DOI archive (e.g. Zenodo). This was not done here: the archive hosts are not
  reachable from this environment.

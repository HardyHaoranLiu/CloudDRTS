# CloudDRTS: preparation work

Code and notes for the four phases of the CloudDRTS preparation plan (October 26 to December 31, 2026). Each note is a short paper in IEEE conference format (`docs/*.tex`, compiled to `docs/*.pdf`).

| Phase | Task | Code | Note |
|---|---|---|---|
| 1 | Kundur two-area system: fault, linearization, inter-area and local modes, comparison with Kundur and with Klein et al. | `src/kundur_system.py`, `src/machines.py`, `src/phase1_*.py` | [docs/phase1_note.pdf](docs/phase1_note.pdf) |
| 2 | Nodal EMT solver (R, L, C, voltage source, switch, lossless Bergeron line), validated at 50 µs and 5 µs | `src/emt.py`, `src/phase2_validation.py` | [docs/phase2_note.pdf](docs/phase2_note.pdf) |
| 3 | Split at the Bergeron line, interface delay, Thevenin equivalents, proposal for the Kundur system | `src/phase3_partition.py` | [docs/phase3_note.pdf](docs/phase3_note.pdf) |
| 4 | Per-step profiling of the solver, questions for the first week | `src/phase4_profile.py` | [docs/phase4_note.pdf](docs/phase4_note.pdf) |

Figures (PDF and PNG) and tables are in `results/`.

## Running

Python 3.12:

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
cd src
../.venv/bin/python run_all.py
```

`run_all.py` runs all scripts in about a minute and rewrites `results/`; each script can also be run on its own from `src/`. numba is only used for one timing in Phase 4, which is skipped if it is not installed. On Windows use `.venv\Scripts\python`.

## Independent check with ANDES (optional)

`src/phase1_andes.py` builds the same systems in [ANDES](https://github.com/CURENT/andes) and compares its eigenvalue results with those of Phase 1. ANDES needs its own environment:

```bash
python -m venv .venv-andes
.venv-andes/bin/pip install -r requirements-andes.txt
cd src
../.venv-andes/bin/python phase1_andes.py      # after phase1_modes.py and phase1_klein.py
```

## Building the notes

The notes use the IEEEtran class and siunitx (both in TeX Live and on Overleaf) and take their figures from `results/`:

```bash
cd docs
pdflatex phase1_note.tex && pdflatex phase1_note.tex
```

and likewise for phases 2 to 4. The Phase 4 numbers are timings and change slightly if `phase4_profile.py` is run again.

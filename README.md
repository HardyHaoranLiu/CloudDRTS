# CloudDRTS: preparation work

Code and notes for the four phases of the CloudDRTS preparation plan (October 26 to December 31, 2026).

| Phase | Task | Code | Note |
|---|---|---|---|
| 1 | Kundur two-area system: fault, linearization, inter-area and local modes, comparison with Kundur and with Klein et al. | `src/kundur_system.py`, `src/machines.py`, `src/phase1_*.py` | [docs/phase1_note.md](docs/phase1_note.md) |
| 2 | Nodal EMT solver (R, L, C, voltage source, switch, lossless Bergeron line), validated at 50 µs and 5 µs | `src/emt.py`, `src/phase2_validation.py` | [docs/phase2_note.md](docs/phase2_note.md) |
| 3 | Split at the Bergeron line, interface delay, Thevenin equivalents, proposal for the Kundur system | `src/phase3_partition.py` | [docs/phase3_note.md](docs/phase3_note.md) |
| 4 | Per-step profiling of the solver, questions for the first week | `src/phase4_profile.py` | [docs/phase4_note.md](docs/phase4_note.md) |

Figures and tables are in `results/`.

## Running

Python 3.12:

```bash
python -m venv .venv
.venv/bin/pip install -r requirements.txt
cd src
../.venv/bin/python run_all.py
```

`run_all.py` runs all scripts in about a minute and rewrites `results/`; each script can also be run on its own from `src/`. numba is only used for one timing in Phase 4, which is skipped if it is not installed. On Windows use `.venv\Scripts\python`.

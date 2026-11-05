"""Run the four onboarding experiments and sensitivity studies in dependency order."""
from pathlib import Path
import os
import subprocess
import sys

def main():
    root = Path(__file__).resolve().parent
    env = dict(os.environ, OPENBLAS_NUM_THREADS="1")
    scripts = (
        "kundur.py", "experiments.py", "waveform.py", "extended.py",
        "grid_cuts.py", "audit_checks.py", "literature_checks.py",
        "event_validation.py", "interface_attribution.py",
        "working_point_validation.py", "detailed_machine.py", "textbook_reproduction.py",
        "phase1_sensitivity.py",
        "detailed_transient.py",
        "dynamic_thevenin.py",
    )
    for index, script in enumerate(scripts, 1):
        print(f"[{index}/{len(scripts)}] {script}", flush=True)
        subprocess.run([sys.executable, str(root / script)], cwd=root, env=env, check=True)
    print(f"Results: {root / 'results'}")

if __name__ == "__main__":
    main()

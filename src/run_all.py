"""Run every script in order. Results go to ../results."""
from pathlib import Path
import subprocess
import sys

SCRIPTS = ["phase1_modes.py", "phase1_textbook.py", "phase1_klein.py",
           "phase2_validation.py", "phase3_partition.py", "phase4_profile.py"]

if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    for script in SCRIPTS:
        print(f"--- {script}", flush=True)
        subprocess.run([sys.executable, script], cwd=here, check=True)

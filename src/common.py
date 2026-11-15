"""Paths, file helpers and plot style shared by all scripts."""
from pathlib import Path
import csv
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RESULTS = Path(__file__).resolve().parents[1] / "results"
RESULTS.mkdir(exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 8.5,
    "axes.labelsize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.1,
    "lines.markersize": 4,
    "legend.frameon": False,
    "axes.prop_cycle": plt.cycler(color=["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9"]),
})


def save_csv(name, header, rows):
    with (RESULTS / f"{name}.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)


def save_json(name, data):
    def plain(x):
        if isinstance(x, (np.floating, np.integer, np.bool_)):
            return x.item()
        if isinstance(x, np.ndarray):
            return x.tolist()
        raise TypeError(type(x))
    (RESULTS / f"{name}.json").write_text(json.dumps(data, indent=2, default=plain))


def panel_labels(axes):
    for k, ax in enumerate(axes):
        ax.text(-0.02, 1.03, f"({chr(97 + k)})", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8.5)


def save_figure(fig, name):
    fig.savefig(RESULTS / f"{name}.png", dpi=200, bbox_inches="tight")
    plt.close(fig)


def rel_rms(x, ref, axis=0):
    """Relative RMS error in percent."""
    return 100 * np.sqrt(np.sum((x - ref) ** 2, axis=axis) / np.sum(ref ** 2, axis=axis))

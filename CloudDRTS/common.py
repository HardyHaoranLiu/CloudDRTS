from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
from cycler import cycler
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Nimbus Sans", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "legend.fontsize": 7,
    "axes.linewidth": .65, "axes.spines.top": True, "axes.spines.right": True,
    "axes.grid": False, "legend.frameon": False, "lines.linewidth": 1.15,
    "lines.markersize": 4, "xtick.direction": "out", "ytick.direction": "out",
    "xtick.major.width": .65, "ytick.major.width": .65,
    "axes.prop_cycle": cycler(color=["#0072B2","#D55E00","#009E73","#AA4499"]),
    "svg.fonttype": "none", "pdf.fonttype": 42
})

def table(name, fields, rows):
    with (OUT / (name + ".csv")).open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerows(rows)

def figure(fig, name):
    for k, ax in enumerate(fig.axes):
        for side in ["left", "right", "top", "bottom"]:
            ax.spines[side].set_visible(True)
        if len(fig.axes) > 1:
            ax.annotate(chr(97+k), (0,1), xycoords="axes fraction", xytext=(-2,5),
                        textcoords="offset points", ha="right", fontweight="bold", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / (name + ".svg"), bbox_inches="tight")
    fig.savefig(OUT / (name + ".pdf"), bbox_inches="tight")
    fig.savefig(OUT / (name + ".png"), dpi=320, bbox_inches="tight")
    plt.close(fig)

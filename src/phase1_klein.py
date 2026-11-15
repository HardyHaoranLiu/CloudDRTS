"""Phase 1: comparison with Klein, Rogers and Kundur (1991).

The paper tabulates the inter-area mode of the symmetric zero-transfer system
(1400 MW generation and 1367 MW constant-impedance load per area):
  Table 2, manual excitation, one tie circuit:        0.481 Hz, zeta 0.023
  Table 1, static exciters, three/two/one circuits:   0.748/0.661/0.513 Hz,
                                                      zeta 0.018/0.011/0.002
Its Fig. 1 shows the units connected directly to the area buses, while the
textbook network has 0.15 p.u. step-up transformers. All four cases are
computed for both, together with a sweep of the transformer reactance and of
the shunt compensation at buses 7 and 9.
"""
import numpy as np

from common import plt, save_csv, save_json, save_figure
from kundur_system import symmetric_case, mode_label
from machines import SixState

PAPER = {("manual", 1): (0.481, 2.3), ("fast", 1): (0.513, 0.2),
         ("fast", 2): (0.661, 1.1), ("fast", 3): (0.748, 1.8)}
NO_TRANSFORMER = 1e-4
TRANSFORMER_X = (NO_TRANSFORMER, 0.025, 0.05, 0.075, 0.1, 0.125, 0.15)
SHUNT = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5)


def inter_area(case, exciter):
    model = SixState(case)
    if exciter == "manual":
        values, vectors = model.modes()
        shapes = vectors[4:8]
    else:
        values, shapes = model.exciter_modes(KA=200.0, TR=0.01)
    labels = [mode_label(shapes[:, j]) for j in range(3)]
    if labels.count("inter-area") != 1:
        raise RuntimeError(f"Unexpected labels {labels}.")
    z = values[labels.index("inter-area")]
    return z.imag / (2 * np.pi), -100 * z.real / abs(z)


def crossing(xs, fs, target):
    for x0, f0, x1, f1 in zip(xs, fs, xs[1:], fs[1:]):
        if (f0 - target) * (f1 - target) <= 0 and f0 != f1:
            return x0 + (target - f0) * (x1 - x0) / (f1 - f0)
    return None


def main():
    rows, cases = [], []
    for (exciter, ties), (f_ref, z_ref) in PAPER.items():
        sweep = [(xt, *inter_area(symmetric_case(xt, ties), exciter)) for xt in TRANSFORMER_X]
        rows += [["transformer", exciter, ties, xt, 200, f, z, f_ref, z_ref] for xt, f, z in sweep]
        xs, fs, zs = zip(*sweep)
        cases.append(dict(exciter=exciter, tie_circuits=ties, paper_Hz=f_ref, paper_damping_pct=z_ref,
                          textbook_network=[fs[-1], zs[-1]], no_transformer=[fs[0], zs[0]],
                          frequency_bracketed=min(fs) <= f_ref <= max(fs),
                          damping_bracketed=min(zs) <= z_ref <= max(zs),
                          transformer_x_at_paper_frequency=crossing(xs, fs, f_ref)))
    shunt = []
    for q in SHUNT:
        f, z = inter_area(symmetric_case(0.15, 1, q), "manual")
        rows.append(["shunt", "manual", 1, 0.15, 100 * q, f, z, *PAPER[("manual", 1)]])
        shunt.append(f)
    save_csv("phase1_klein", ["sweep", "exciter", "tie_circuits", "transformer_x_pu",
                              "shunt_Mvar_per_area", "inter_area_Hz", "inter_area_damping_pct",
                              "paper_Hz", "paper_damping_pct"], rows)
    save_json("phase1_klein", dict(cases=cases, shunt_sweep_range_Hz=[min(shunt), max(shunt)]))

    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    names = ["Manual, 1 tie", "Static exciter, 1 tie", "Static exciter, 2 ties", "Static exciter, 3 ties"]
    for k, (c, name) in enumerate(zip(cases, names)):
        r = [x for x in rows if x[0] == "transformer" and x[1] == c["exciter"] and x[2] == c["tie_circuits"]]
        ax.plot([x[3] for x in r], [x[5] for x in r], "o-", color=f"C{k}", mfc="white", ms=3.5, label=name)
        ax.axhline(c["paper_Hz"], color=f"C{k}", lw=0.7, ls=":")
        if c["transformer_x_at_paper_frequency"] is not None:
            ax.plot(c["transformer_x_at_paper_frequency"], c["paper_Hz"], "D", color=f"C{k}", ms=4)
    ax.axvline(0.15, color="0.5", lw=0.7, ls="--")
    ax.text(0.148, 0.795, "textbook network", ha="right", va="top", fontsize=6.5, color="0.35")
    ax.set(xlabel="Step-up transformer reactance (p.u., 900 MVA base)",
           ylabel="Inter-area mode (Hz)", xlim=(-0.005, 0.16), ylim=(0.39, 0.80))
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.465), fontsize=6.5)
    save_figure(fig, "phase1_klein")
    for c in cases:
        print(c)
    print("shunt sweep", min(shunt), max(shunt))


if __name__ == "__main__":
    main()

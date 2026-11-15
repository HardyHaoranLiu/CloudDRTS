"""Phase 1: comparison of the inter-area mode with Klein, Rogers and Kundur (1991).

The paper tabulates the inter-area mode of its symmetric zero-transfer system
(1400 MW generation and 1367 MW constant-impedance load per area, all H = 6.5 s):
  Table 2, manual excitation, one tie circuit:        0.481 Hz, zeta 0.023
  Table 1, static exciters, three/two/one circuits:   0.748/0.661/0.513 Hz,
                                                      zeta 0.018/0.011/0.002
Its Fig. 1 connects the units directly to the area buses (no step-up
transformers) with the line data of the textbook network. The shunt
capacitors are only described as added as necessary; 200 Mvar per area is
used here. The four cases are computed for this network, together with the
sensitivity to data the paper does not give, and the textbook network
(with 0.15 p.u. step-up transformers) for reference.
"""
import numpy as np

from common import save_csv, save_json
from kundur_system import symmetric_case, mode_label
from machines import SixState

PAPER = {("manual", 1): (0.481, 2.3), ("fast", 1): (0.513, 0.2),
         ("fast", 2): (0.661, 1.1), ("fast", 3): (0.748, 1.8)}
DIRECT = 1e-4                      # transformer reactance standing for a direct connection
VARIANTS = {"paper network": {},
            "shunt 100 Mvar per area": dict(shunt=1.0),
            "shunt 300 Mvar per area": dict(shunt=3.0),
            "no line charging": dict(b_km=1e-9),
            "no line resistance": dict(z_km=1e-3j)}


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
    j = labels.index("inter-area")
    z = values[j]
    shape = shapes[:, j] / shapes[np.argmax(abs(shapes[:, j])), j]
    return z.imag / (2 * np.pi), -100 * z.real / abs(z), shape


def main():
    rows, cases = [], []
    for name, kw in VARIANTS.items():
        for (exciter, ties), (f_ref, z_ref) in PAPER.items():
            f, z, _ = inter_area(symmetric_case(DIRECT, ties, **kw), exciter)
            rows.append([name, exciter, ties, f, z, f_ref, z_ref, 100 * (f / f_ref - 1)])
    for (exciter, ties), (f_ref, z_ref) in PAPER.items():
        f, z, _ = inter_area(symmetric_case(0.15, ties), exciter)
        rows.append(["textbook network", exciter, ties, f, z, f_ref, z_ref, 100 * (f / f_ref - 1)])
    save_csv("phase1_klein", ["network", "exciter", "tie_circuits", "inter_area_Hz",
                              "inter_area_damping_pct", "paper_Hz", "paper_damping_pct",
                              "frequency_difference_pct"], rows)

    # Mode shape of the paper network, manual excitation, one tie (cf. Fig. 2 of the paper).
    _, _, shape = inter_area(symmetric_case(DIRECT, 1), "manual")
    base = [r for r in rows if r[0] == "paper network"]
    save_json("phase1_klein", dict(
        paper_network_frequency_difference_pct=[min(r[7] for r in base), max(r[7] for r in base)],
        paper_network_max_damping_difference_points=max(abs(r[4] - r[6]) for r in base),
        sensitivity_frequency_change_pct={name: max(abs(r[3] / b[3] - 1) * 100 for r, b in
                                                    zip([x for x in rows if x[0] == name], base))
                                          for name in VARIANTS if name != "paper network"},
        mode_shape_G1_to_G4=dict(amplitude=abs(shape).tolist(), phase_deg=np.angle(shape, deg=True).tolist())))
    for r in rows:
        print(f"{r[0]:24s} {r[1]:6s} {r[2]}  {r[3]:.3f} Hz {r[4]:5.2f} %   paper {r[5]} / {r[6]}   {r[7]:+.1f} %")
    print("shape", np.round(abs(shape), 3), np.round(np.angle(shape, deg=True), 1))


if __name__ == "__main__":
    main()

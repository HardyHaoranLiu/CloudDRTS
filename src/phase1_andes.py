"""Phase 1: independent check of the modal results with ANDES.

The same systems are built in ANDES (GENROU machines without saturation,
same network and load models) and analysed with its eigenvalue routine:
  textbook   Kundur Example 12.6, manual excitation, active load constant
             current and reactive load constant impedance
  klein      network of Klein et al. (no step-up transformers), symmetric
             zero-transfer case, constant-impedance loads, H = 6.5 s,
             manual excitation (one tie) and static exciters (1-3 ties)
The ANDES results are compared with results/phase1_modes.csv and
results/phase1_klein.csv. Requires the packages in requirements-andes.txt
(run phase1_modes.py and phase1_klein.py first).
"""
import csv
import json
from pathlib import Path

import numpy as np
import andes

RESULTS = Path(__file__).resolve().parents[1] / "results"
KM = (1e-4, 1e-3, 1.75e-3)
LINES = [(5, 6, 25), (6, 7, 10), (7, 8, 110), (8, 9, 110), (9, 10, 10), (10, 11, 25)]
UNITS = [(1, 5), (2, 6), (3, 11), (4, 10)]


def build(transformer_x, ties, shunt, loads, pv, H, p2i, exciter):
    andes.config_logger(stream_level=50)
    ss = andes.System(no_output=True, default_config=True)
    for i in range(1, 12):
        ss.add("Bus", dict(idx=i, name=f"B{i}", Vn=20.0 if i <= 4 else 230.0))
    n = 0
    for a, b, km in LINES:
        for _ in range(ties if km == 110 else 1):
            n += 1
            ss.add("Line", dict(idx=f"L{n}", bus1=a, bus2=b, r=KM[0] * km, x=KM[1] * km,
                                b=KM[2] * km, Vn1=230.0, Vn2=230.0))
    for k, (g, hv) in enumerate(UNITS):
        ss.add("Line", dict(idx=f"T{k+1}", bus1=g, bus2=hv, r=0.0, x=transformer_x / 9, b=0.0,
                            trans=1, tap=1.0, Vn1=20.0, Vn2=230.0))
    for bus, b in zip((7, 9), shunt):
        ss.add("Shunt", dict(idx=f"C{bus}", bus=bus, Vn=230.0, Sn=100.0, b=b, g=0.0))
    for bus, (p, q) in zip((7, 9), loads):
        ss.add("PQ", dict(idx=f"P{bus}", bus=bus, Vn=230.0, p0=p, q0=q))
    ss.add("Slack", dict(idx="G1", bus=1, Vn=20.0, v0=1.03, a0=0.0, p0=7.0))
    for bus, p, v in pv:
        ss.add("PV", dict(idx=f"G{bus}", bus=bus, Vn=20.0, p0=p, v0=v))
    for k in range(4):
        ss.add("GENROU", dict(idx=f"M{k+1}", bus=k + 1, gen=f"G{k+1}", Sn=900.0, Vn=20.0, fn=60.0,
                              M=2 * H[k], D=0.0, ra=0.0025, xl=0.2, xd=1.8, xq=1.7, xd1=0.3,
                              xq1=0.55, xd2=0.25, xq2=0.25, Td10=8.0, Td20=0.03, Tq10=0.4,
                              Tq20=0.05, S10=0.0, S12=0.0))
        if exciter:
            ss.add("EXST1", dict(idx=f"E{k+1}", syn=f"M{k+1}", TR=0.01, VIMAX=99.0, VIMIN=-99.0,
                                 TC=0.0, TB=0.0, KA=200.0, TA=0.0, VRMAX=99.0, VRMIN=-99.0,
                                 KC=0.0, KF=0.0, TF=1.0))
    ss.PQ.config.p2p, ss.PQ.config.p2i, ss.PQ.config.p2z = 0.0, float(p2i), float(not p2i)
    ss.PQ.config.q2q, ss.PQ.config.q2i, ss.PQ.config.q2z = 0.0, 0.0, 1.0
    ss.setup()
    return ss


def rotor_modes(ss):
    ss.PFlow.run()
    ss.EIG.run()
    mu = np.asarray(ss.EIG.mu)
    _, _, right, left = ss.EIG.calc_pfactor()
    pf = np.abs(left) * np.abs(right)                 # states x modes
    rotor = np.r_[ss.GENROU.delta.a, ss.GENROU.omega.a]
    share = pf[rotor].sum(0) / pf.sum(0)
    ix = sorted([i for i in range(len(mu)) if 1 < mu[i].imag < 12 and share[i] > 0.3],
                key=lambda i: mu[i].imag)
    if len(ix) != 3:
        raise RuntimeError(f"Expected three rotor modes, found {len(ix)}.")
    return mu[ix], mu


def symmetric(transformer_x, ties, exciter):
    pv = [(2, 7.0, 1.01), (3, 7.0, 1.03), (4, 7.0, 1.01)]
    for _ in range(30):
        ss = build(transformer_x, ties, (2.0, 2.0), ((13.67, 1.0), (13.67, 1.0)), pv,
                   (6.5,) * 4, p2i=False, exciter=exciter)
        ss.PFlow.run()
        slack = float(ss.Slack.p.v[0])
        if abs(slack - pv[1][1]) < 1e-10:
            break
        pv[1] = (3, slack, 1.03)              # same output for G1 and G3: zero transfer
    return ss


def fz(z):
    return z.imag / (2 * np.pi), -100 * z.real / abs(z)


def main():
    rows = []
    ours = {r["mode"]: (float(r["frequency_Hz"]), float(r["damping_pct"]))
            for r in csv.DictReader((RESULTS / "phase1_modes.csv").open())}
    ss = build(0.15, 2, (2.0, 3.5), ((9.67, 1.0), (17.67, 1.0)),
               [(2, 7.0, 1.01), (3, 7.19, 1.03), (4, 7.0, 1.01)],
               (6.5, 6.5, 6.175, 6.175), p2i=True, exciter=False)
    modes, spectrum = rotor_modes(ss)
    for name, z in zip(["inter-area", "area-1", "area-2"], modes):
        f, d = fz(z)
        rows.append(["textbook", name, f, d, *ours[name], f - ours[name][0], d - ours[name][1]])
    slow = sorted(float(z.real) for z in spectrum if abs(z.imag) < 1e-9 and z.real > 1e-3)

    klein = {(r["exciter"], int(r["tie_circuits"])): (float(r["inter_area_Hz"]), float(r["inter_area_damping_pct"]))
             for r in csv.DictReader((RESULTS / "phase1_klein.csv").open()) if r["network"] == "paper network"}
    paper = {("manual", 1): 0.481, ("fast", 1): 0.513, ("fast", 2): 0.661, ("fast", 3): 0.748}
    for (exc, ties), ref in paper.items():
        modes, _ = rotor_modes(symmetric(1e-4, ties, exc == "fast"))
        f, d = fz(modes[0])
        o = klein[(exc, ties)]
        rows.append([f"klein, {exc}, {ties} tie(s)", "inter-area", f, d, *o, f - o[0], d - o[1]])
        rows[-1].append(100 * (f / ref - 1))

    with (RESULTS / "phase1_andes.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["case", "mode", "andes_Hz", "andes_damping_pct", "this_model_Hz",
                    "this_model_damping_pct", "diff_Hz", "diff_damping_points", "andes_vs_paper_pct"])
        w.writerows(rows)
    summary = dict(andes_version=andes.__version__,
                   max_frequency_difference_Hz=max(abs(r[6]) for r in rows),
                   max_damping_difference_points=max(abs(r[7]) for r in rows),
                   andes_vs_klein_frequency_pct=[r[8] for r in rows if len(r) > 8],
                   textbook_positive_real_eigenvalues=slow)
    (RESULTS / "phase1_andes.json").write_text(json.dumps(summary, indent=2))
    for r in rows:
        print(r)
    print(summary)


if __name__ == "__main__":
    main()

"""Phase 1: comparison with Kundur's Example 12.6 (Tables E12.3 and E12.4).

The example uses the equivalent-circuit machine model with saturation
(Asat 0.015, Bsat 9.6, psiT1 0.9), linearized with the incremental saturation
factor. machines.EquivalentCircuit implements that model; without saturation
it has the same eigenvalues as machines.SixState.
"""
import numpy as np

from common import save_csv, save_json
from kundur_system import textbook_case, mode_label
from machines import SixState, EquivalentCircuit

# Table E12.3 (manual excitation), one entry per real eigenvalue or complex pair,
# without the pair at the origin.
TABLE_E123 = [-0.096, -0.111 + 3.43j, -0.117, -0.265, -0.276, -0.492 + 6.82j, -0.506 + 7.02j,
              -3.428, -4.139, -5.287, -5.303, -31.03, -32.45, -34.07, -35.53,
              -37.89 + 0.142j, -38.01 + 0.038j]
# Table E12.4, rotor modes (inter-area, area 1, area 2).
TABLE_E124 = {"thyristor": (0.031 + 3.84j, -0.490 + 7.15j, -0.496 + 7.35j),
              "tgr": (0.123 + 3.46j, -0.450 + 6.86j, -0.462 + 7.05j),
              "pss": (-0.501 + 3.77j, -1.826 + 8.05j, -1.895 + 8.35j)}
# Frequencies (Hz) and damping ratios (%) printed with Table E12.3.
BOOK = {"inter-area": (0.545, 3.2), "area-1": (1.087, 7.2), "area-2": (1.117, 7.2)}
MODES = ["inter-area", "area-1", "area-2"]


def f_zeta(z):
    return z.imag / (2 * np.pi), -100 * z.real / abs(z)


def main():
    case = textbook_case()
    six = SixState(case)

    # Same eigenvalues without saturation?
    a = np.sort_complex(EquivalentCircuit(case, "none").spectrum())
    b = np.sort_complex(six.eigen(1e-5)[0])
    equivalence = float(np.max(abs(a - b)))

    book = EquivalentCircuit(case, "incremental")
    spectrum = book.spectrum()
    rest = [z for z in spectrum if abs(z) > 1e-3 and z.imag >= 0]
    rows, e123 = [], []
    for ref in TABLE_E123:
        ours = min(rest, key=lambda z: abs(z - ref))
        rest.remove(ours)
        e123.append(abs(ours - ref))
        rows.append(["E12.3", "manual", ref.real, ref.imag, ours.real, ours.imag, abs(ours - ref)])
    e124 = {}
    for kind, refs in TABLE_E124.items():
        modes = book.rotor_modes(kind)
        diffs = []
        for label, ref in zip(MODES, refs):
            diffs.append(abs(modes[label] - ref))
            rows.append(["E12.4", f"{kind}/{label}", ref.real, ref.imag,
                         modes[label].real, modes[label].imag, diffs[-1]])
        e124[kind] = max(diffs)
    save_csv("phase1_textbook_eigenvalues",
             ["table", "case", "book_real", "book_imag", "real", "imag", "abs_difference"], rows)

    # Mode table: this model without saturation, the book's method, the book.
    values, vectors = six.modes()
    six_modes = {mode_label(v[4:8]): z for z, v in zip(values, vectors.T)}
    manual = book.rotor_modes("manual")
    table = []
    for label in MODES:
        f0, z0 = f_zeta(six_modes[label])
        f1, z1 = f_zeta(manual[label])
        table.append([label, f0, z0, f1, z1, *BOOK[label]])
    save_csv("phase1_mode_table", ["mode", "f_no_sat_Hz", "zeta_no_sat_pct", "f_book_method_Hz",
                                   "zeta_book_method_pct", "f_book_Hz", "zeta_book_pct"], table)
    thyristor = book.rotor_modes("thyristor")["inter-area"]
    save_json("phase1_textbook", dict(
        six_state_vs_equivalent_circuit_max_eigenvalue_difference=equivalence,
        table_E123_max_abs_difference=max(e123),
        table_E123_entries=len(e123),
        table_E124_max_abs_difference=e124,
        thyristor_inter_area_eigenvalue=[thyristor.real, thyristor.imag],
        max_real_eigenvalue_book_method=float(max(z.real for z in spectrum if abs(z) > 1e-3))))
    for r in table:
        print("{:>10}: {:.4f} Hz {:.2f} % | {:.4f} Hz {:.2f} % | {} Hz {} %".format(*r))
    print("E12.3 max diff", max(e123), "E12.4", e124, "equivalence", equivalence)


if __name__ == "__main__":
    main()

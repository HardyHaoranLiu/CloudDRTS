"""Reproduce the inter-area results of Klein, Rogers and Kundur [2] (Tables 1-2).

[2] reports, for the symmetric zero-transfer system (1400 MW generation and
1367 MW constant-impedance load per area, detailed machines of its Appendix A):
  Table 2, manual excitation, one tie circuit:     0.481 Hz, damping ratio 0.023
  Table 1, fast static exciters, 3/2/1 tie circuits: 0.748/0.661/0.513 Hz,
                                                    0.018/0.011/0.002
[2] gives no local-mode values. Its Fig. 1 draws the units directly on the
area buses and gives no step-up transformer data, while the textbook network
[1] has 0.15 p.u. (900 MVA base) step-up transformers. Both representations
are computed, plus a sweep of that reactance and of the shunt compensation.
This is a declared sweep; nothing is fitted or adopted from it.
"""
import json
from common import np, table, OUT
import kundur as kd
from literature_checks import symmetric
from detailed_machine import Detailed, static_exciter_modes

# [2] Table 1 (fast static exciter, no TGR) and Table 2 (manual): (Hz, damping %).
KLEIN = {('manual', 1): (.481, 2.3), ('fast', 1): (.513, .2),
         ('fast', 2): (.661, 1.1), ('fast', 3): (.748, 1.8)}
KA, TR = 200., .01           # [2] Appendix A fast static exciter; TC = TB = 0 (no TGR)
NO_TRANSFORMER = 1e-4        # p.u. on the machine base; stands for "no step-up transformer"
TRANSFORMER_X = (NO_TRANSFORMER, .025, .05, .075, .1, .125, .15)
COMPENSATION_PU = (1.0, 1.5, 2.0, 2.5, 3.0, 3.5)

def inter_area(case, exciter, model='six_state'):
    """Inter-area eigenvalue (Hz, damping %) with manual or fast static excitation."""
    if model == 'two_axis':
        value = case.modes()[0][0]
        return value.imag/(2*np.pi), -100*value.real/abs(value)
    m = Detailed(case, 'none')
    if exciter == 'manual':
        values, vectors = m.modes()
        assert kd.mode_label(vectors[4:8, 0]) == 'inter-area'
        value = values[0]
        return value.imag/(2*np.pi), -100*value.real/abs(value)
    values, shapes = static_exciter_modes(m, KA, TR)
    labels = [kd.mode_label(shapes[:, j]) for j in range(3)]
    if labels.count('inter-area') != 1:
        raise RuntimeError(f'Unexpected labels {labels}.')
    value = values[labels.index('inter-area')]
    return value.imag/(2*np.pi), -100*value.real/abs(value)

def crossing(xs, fs, target):
    """X_T at which the (monotone) frequency curve reaches the [2] value, else None."""
    for (x0, f0), (x1, f1) in zip(zip(xs, fs), zip(xs[1:], fs[1:])):
        if (f0-target)*(f1-target) <= 0 and f0 != f1:
            return x0+(target-f0)*(x1-x0)/(f1-f0)
    return None

def main():
    rows = []
    for (exciter, ties), (ref_f, ref_z) in KLEIN.items():
        for xt in TRANSFORMER_X:
            f, z = inter_area(symmetric(2.0, ties, xt), exciter)
            rows.append(['transformer_x', exciter, ties, xt, 200, 'six_state', f, z, ref_f, ref_z])
    for qc in COMPENSATION_PU:
        for model in ('six_state', 'two_axis'):
            f, z = inter_area(symmetric(qc, 1, .15), 'manual', model)
            rows.append(['compensation', 'manual', 1, .15, qc*100, model, f, z, *KLEIN[('manual', 1)]])
    table('phase1_sensitivity', ['sweep', 'exciter', 'tie_circuits', 'transformer_x_pu_machine',
          'compensation_Mvar_per_area', 'machine_model', 'inter_area_Hz', 'inter_area_damping_pct',
          'ref2_Hz', 'ref2_damping_pct'], rows)

    cases = []
    for (exciter, ties), (ref_f, ref_z) in KLEIN.items():
        rr = [r for r in rows if r[0] == 'transformer_x' and r[1] == exciter and r[2] == ties]
        xs = [r[3] for r in rr]; fs = [r[6] for r in rr]; zs = [r[7] for r in rr]
        textbook = rr[-1]; none = rr[0]
        cases.append(dict(exciter=exciter, tie_circuits=ties, ref2_Hz=ref_f, ref2_damping_pct=ref_z,
            textbook_network_Hz=textbook[6], textbook_network_damping_pct=textbook[7],
            no_transformer_Hz=none[6], no_transformer_damping_pct=none[7],
            frequency_bracketed=bool(min(fs) <= ref_f <= max(fs)),
            damping_bracketed=bool(min(zs) <= ref_z <= max(zs)),
            frequency_crossing_transformer_x_pu=crossing(xs, fs, ref_f)))
    comp = [r for r in rows if r[0] == 'compensation' and r[5] == 'six_state']
    two = [r for r in rows if r[0] == 'compensation' and r[5] == 'two_axis' and r[4] == 200][0]
    base = [r for r in comp if r[4] == 200][0]
    crossings = [c['frequency_crossing_transformer_x_pu'] for c in cases]
    summary = dict(
        source='Klein, Rogers, Kundur 1991, Tables 1 and 2; values checked against the paper.',
        cases=cases,
        all_frequencies_bracketed=all(c['frequency_bracketed'] for c in cases),
        all_dampings_bracketed=all(c['damping_bracketed'] for c in cases),
        frequency_crossing_range_transformer_x_pu=[min(crossings), max(crossings)] if None not in crossings else None,
        v7_case=dict(six_state_Hz=base[6], six_state_damping_pct=base[7], two_axis_Hz=two[6]),
        compensation_range_Hz=[min(r[6] for r in comp), max(r[6] for r in comp)],
        machine_model_effect_at_v7_case_Hz=abs(base[6]-two[6]),
        note='Sweep only; no parameter is fitted or adopted.')
    (OUT/'phase1_sensitivity.json').write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()

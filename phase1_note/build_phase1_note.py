"""Build the Phase 1 two-page notes (English and Chinese) in IEEE conference layout.

Every number printed in the notes comes from a result file through
build/values.tex; statements in the text that depend on results are checked
by the assertions in collect(), so a changed result stops the build instead of
leaving stale prose.

Inputs (create them first):
  CloudDRTS/run_all.py                    -> CloudDRTS/results/
  external_run/run_reproduction.py        -> external_run/results/
Toolchain: XeLaTeX with IEEEtran, xeCJK (ctex) and siunitx; fonts TeX Gyre
Termes and Noto Serif CJK SC. Python: numpy, matplotlib (pgf backend), PyMuPDF.

Outputs: output/pdf/Phase1_Note_EN.pdf, output/pdf/Phase1_Note_ZH.pdf and
page previews in qa/Phase1_Note_{EN,ZH}/.
"""
from pathlib import Path
import csv, hashlib, json, shutil, subprocess
import numpy as np
import matplotlib
matplotlib.use('pgf')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
RES = ROOT/'CloudDRTS'/'results'
EXT = ROOT/'external_run'/'results'
BUILD = HERE/'build'
PDF_OUT = ROOT/'output'/'pdf'
H = np.array([6.5, 6.5, 6.175, 6.175])
FAULT = (1.0, 1.08)

PREAMBLE = (r'\usepackage{mathptmx}\usepackage[no-math]{fontspec}'
            r'\setmainfont{TeX Gyre Termes}\usepackage{xeCJK}'
            r'\setCJKmainfont{Noto Serif CJK SC}')
plt.rcParams.update({
    'pgf.texsystem': 'xelatex', 'pgf.rcfonts': False, 'pgf.preamble': PREAMBLE,
    'font.family': 'serif', 'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8,
    'legend.fontsize': 7, 'xtick.labelsize': 7, 'ytick.labelsize': 7,
    'axes.linewidth': .5, 'lines.linewidth': .9, 'xtick.major.width': .5,
    'ytick.major.width': .5, 'xtick.direction': 'in', 'ytick.direction': 'in',
    'legend.frameon': False, 'savefig.pad_inches': .02,
    'axes.prop_cycle': plt.cycler(color=['#0072B2', '#D55E00', '#009E73', '#CC79A7']),
})
COLUMN = 3.45  # inches; IEEEtran conference \columnwidth is 3.5 in

TEXT = {
    'en': dict(modes=['Inter-area', 'Local, area 1', 'Local, area 2'],
               angle='Rotor angle\nvs. COI (deg)', tie='Tie power\n7→8 (MW)',
               coi='COI frequency\ndeviation (Hz)', time='Time (s)',
               ours='This note', andes='ANDES 1.9.3',
               xt='Step-up transformer reactance (p.u., 900 MVA base)', freq='Inter-area mode (Hz)',
               cases=['Manual, 1 tie', 'Fast, 1 tie', 'Fast, 2 ties', 'Fast, 3 ties'],
               ref='[2]', textbook='[1] network'),
    'zh': dict(modes=['区间模态', '区域 1 本地模态', '区域 2 本地模态'],
               angle='相对惯量中心\n功角 (°)', tie='联络线功率\n7→8 (MW)',
               coi='惯量中心\n频率偏差 (Hz)', time='时间 (s)',
               ours='本文', andes='ANDES 1.9.3',
               xt='升压变压器电抗 (p.u.，900 MVA 基准)', freq='区间模态频率 (Hz)',
               cases=['手动励磁，1 回', '快速励磁，1 回', '快速励磁，2 回', '快速励磁，3 回'],
               ref='[2]', textbook='[1] 网络'),
}

def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

def need(path):
    if not path.exists():
        raise SystemExit(f'Missing {path.relative_to(ROOT)}. Run CloudDRTS/run_all.py and '
                         'external_run/run_reproduction.py first (see README).')
    return path

def num(x, digits=3):
    return rf'\num{{{x:.{digits}f}}}'

def sci(x, digits=2):
    return rf'\num{{{x:.{digits-1}e}}}'

def collect():
    modes = [r for r in rows(need(RES/'detailed_modes.csv')) if r['case'] == 'textbook']
    base = [r for r in modes if r['saturation'] == 'none']
    sat = [r for r in modes if r['saturation'] == 'kundur_exponential']
    assert [r['mode_label'] for r in sat] == ['inter-area', 'area-1', 'area-2']
    book = rows(need(RES/'textbook_eigenvalue_comparison.csv'))
    thy = [r for r in book if r['saturation'] == 'kundur_exponential' and r['excitation'] == 'thyristor']
    assert [r['mode_label'] for r in base] == ['inter-area', 'area-1', 'area-2'], 'mode labels changed'
    paper = next(r for r in rows(RES/'detailed_modes.csv')
                 if r['case'] == 'paper_symmetric' and r['saturation'] == 'none' and r['mode_label'] == 'inter-area')
    two_axis = rows(need(RES/'modes.csv'))
    checks = json.loads(need(RES/'detailed_checks.json').read_text())
    textbook_check = next(c for c in checks if c['case'] == 'textbook' and c['saturation'] == 'none')
    check = {(c['case'], c['saturation']): c for c in checks}
    transient = json.loads(need(RES/'detailed_transient_checks.json').read_text())
    sens = json.loads(need(RES/'phase1_sensitivity.json').read_text())
    pf = json.loads(need(RES/'kundur_checks.json').read_text())
    modal = json.loads(need(EXT/'modal_comparison_checks.json').read_text())
    summary = json.loads(need(EXT/'verification_summary.json').read_text())
    fault_csv = need(RES/'detailed_fault.csv')
    assert summary['local_fault_csv_sha256'] == hashlib.sha256(fault_csv.read_bytes()).hexdigest(), \
        'ANDES comparison was made against a different detailed_fault.csv; rerun external_run.'
    assert modal['acceptance']['passed']
    fine = {m['quantity']: m['max_abs'] for m in summary['waveform_metrics'] if m['run'] == 'fine'}

    # Claims made in the text.
    sat_shift = max(abs(float(r['frequency_Hz'])-float(b['frequency_Hz']))
                    for r in modes for b in base if r['mode'] == b['mode'])
    assert sat_shift < .01, 'text says saturation shifts each mode by less than 0.01 Hz'
    assert transient['max_modal_state_error_ratio'] < .005
    assert sens['all_frequencies_bracketed'] and sens['all_dampings_bracketed'], \
        'text says the two network representations bracket every [2] frequency and damping'
    assert sens['compensation_range_Hz'][1] < .481, 'text says compensation alone cannot reach 0.481 Hz'
    assert sens['frequency_crossing_range_transformer_x_pu'] is not None
    assert sens['machine_model_effect_at_v7_case_Hz'] < .002
    assert transient['final_COI_df_Hz'] > .1, 'text says the COI frequency does not return'
    df_tb = [float(r['frequency_error_pct']) for r in base]
    df_sat = [float(r['frequency_error_pct']) for r in sat]
    df_thy = [float(r['frequency_difference_pct']) for r in thy]
    assert max(map(abs, df_tb + df_sat + df_thy)) < 1.5, 'text says frequencies agree with [1] to about 1 %'
    assert all(float(r['sigma_per_s']) > 0 and float(r['ref1_sigma_per_s']) > 0
               for r in thy if r['mode_label'] == 'inter-area'), 'text says both find an unstable inter-area mode'
    excess = [float(r['damping_pct'])-float(r['reference_damping_pct']) for r in sat[1:]]
    assert all(e > 0 for e in excess), 'text says local damping is higher than in [1]'
    pos = check[('textbook', 'none')]['positive_real_eigenvalues']
    pos_andes = modal['positive_real_eigenvalues']
    assert len(pos) == 1 and len(pos_andes) == 1 and abs(pos[0]-pos_andes[0]) < 1e-4, \
        'text says the unsaturated model has one slow positive real eigenvalue, confirmed by ANDES'
    assert not check[('textbook', 'kundur_exponential')]['positive_real_eigenvalues']
    pos_sym = check[('paper_symmetric', 'none')]['positive_real_eigenvalues']
    assert len(pos_sym) == 1 and not check[('paper_symmetric', 'kundur_exponential')]['positive_real_eigenvalues']
    repro = json.loads(need(RES/'textbook_reproduction.json').read_text())
    assert repro['table_E123_entries'] == 17 and repro['table_E123_max_abs_difference'] < .015, \
        'text says all eigenvalues of Table E12.3 are reproduced'
    thy_tgr = max(repro['table_E124_max_abs_difference'][k] for k in ('thyristor', 'tgr'))
    assert thy_tgr < .01 and repro['genrou_vs_circuit_no_saturation_max_eigenvalue_difference'] < 2e-3
    sym = {r['saturation']: float(r['frequency_Hz']) for r in rows(RES/'detailed_modes.csv')
           if r['case'] == 'paper_symmetric' and r['mode_label'] == 'inter-area'}
    symsat = abs(sym['none']-sym['kundur_exponential'])
    assert symsat < .001

    fault = np.genfromtxt(fault_csv, delimiter=',', names=True)
    commit = subprocess.run(['git', 'rev-parse', '--short', 'HEAD'], cwd=ROOT,
                            capture_output=True, text=True).stdout.strip() or 'unknown'
    v = {}
    for key, r in zip(['IA', 'LA', 'LB'], base):
        v['f'+key] = num(float(r['frequency_Hz']), 4)
        v['s'+key] = num(float(r['sigma_per_s']), 3)
        v['z'+key] = num(float(r['damping_pct']), 2)
        v['fref'+key] = num(float(r['reference_Hz']), 3)
        v['zref'+key] = num(float(r['reference_damping_pct']), 1)
        v['df'+key] = rf'\num{{{float(r["frequency_error_pct"]):+.2f}}}'
    for key, r, t in zip(['IA', 'LA', 'LB'], sat, thy):
        v['fs'+key] = num(float(r['frequency_Hz']), 4)
        v['zs'+key] = num(float(r['damping_pct']), 2)
        v['ft'+key] = num(float(t['frequency_Hz']), 2)
        v['ftref'+key] = num(float(t['ref1_frequency_Hz']), 2)
    for key, label in zip(['IA', 'LA', 'LB'], ['inter-area', 'area-1', 'area-2']):
        mode = repro['manual_modes'][label]
        v['fr'+key] = num(mode['frequency_Hz'], 4); v['zr'+key] = num(mode['damping_pct'], 2)
        r = next(x for x in two_axis if x['mode'] == label)
        v['ftwo'+key] = num(float(r['freq_Hz']), 4)
    v.update(
        maxdf=num(max(map(abs, df_tb)), 2), satshift=num(sat_shift, 4),
        maxdfsat=num(max(map(abs, df_sat)), 1), maxdfthy=num(max(map(abs, df_thy)), 1),
        dzlo=num(min(excess), 1), dzhi=num(max(excess), 1),
        posnone=num(pos[0], 4), posandes=num(pos_andes[0], 4), tau=rf'\num{{{round(1/pos[0], -1):.0f}}}',
        possym=num(pos_sym[0], 3), symsat=num(symsat, 4),
        eonetwothree=num(repro['table_E123_max_abs_difference'], 3), thydiff=num(thy_tgr, 3),
        pssdiff=num(repro['table_E124_max_abs_difference']['pss'], 2),
        equiv=sci(repro['genrou_vs_circuit_no_saturation_max_eigenvalue_difference'], 1),
        pfres=sci(pf['pf_residual']), tie=num(fault['P7to8_MW'][0], 1),
        eqres=sci(textbook_check['equilibrium_max_residual']),
        jacdiff=sci(textbook_check['jacobian_halving_max_eigenvalue_change']),
        modalerr=num(100*transient['max_modal_state_error_ratio'], 2),
        prefault=sci(transient['prefault_max_speed_deviation_Hz']),
        vmin=num(transient['minimum_fault_V8_pu'], 4), coidf=num(transient['final_COI_df_Hz'], 2),
        andesdf=sci(modal['matched_frequency_max_abs_difference_Hz']),
        andeseig=sci(modal['matched_eigenvalue_max_abs_difference_per_s']),
        andesv=sci(modal['matched_voltage_max_abs_difference_pu']),
        andesang=sci(fine['COI_angle']), andesfreq=sci(fine['speed_deviation']),
        andesvolt=sci(fine['bus8_voltage']), andestie=sci(fine['tie_power']),
        fsym=num(float(paper['frequency_Hz']), 4), zsym=num(float(paper['damping_pct']), 2),
        dfsym=rf'\num{{{float(paper["frequency_error_pct"]):+.1f}}}',
        frefsym=num(float(paper['reference_Hz']), 3), zrefsym=num(float(paper['reference_damping_pct']), 1),
        fsymtwo=num(sens['v7_case']['two_axis_Hz'], 4),
        complo=num(sens['compensation_range_Hz'][0], 3), comphi=num(sens['compensation_range_Hz'][1], 3),
        xtlo=num(sens['frequency_crossing_range_transformer_x_pu'][0], 3),
        xthi=num(sens['frequency_crossing_range_transformer_x_pu'][1], 3),
        xtfastlo=num(min(c['frequency_crossing_transformer_x_pu'] for c in sens['cases'] if c['exciter'] == 'fast'), 3),
        xtfasthi=num(max(c['frequency_crossing_transformer_x_pu'] for c in sens['cases'] if c['exciter'] == 'fast'), 3),
        modeleffect=num(sens['machine_model_effect_at_v7_case_Hz'], 4),
        commit=commit)
    for c, key in zip(sens['cases'], ['M', 'Fa', 'Fb', 'Fc']):
        v['ref'+key] = rf"{c['ref2_Hz']:.3f}\,/\,{c['ref2_damping_pct']:.1f}"
        v['tb'+key] = rf"{c['textbook_network_Hz']:.3f}\,/\,\num{{{c['textbook_network_damping_pct']:.1f}}}"
        v['nt'+key] = rf"{c['no_transformer_Hz']:.3f}\,/\,\num{{{c['no_transformer_damping_pct']:.1f}}}"
        v['xt'+key] = num(c['frequency_crossing_transformer_x_pu'], 3)
    return dict(values=v, base=base, fault=fault, sens=rows(RES/'phase1_sensitivity.csv'), cases=sens['cases'],
                andes=np.genfromtxt(need(EXT/'matched_fault_fine_waveform.csv'), delimiter=',', names=True))

def write_values(v):
    lines = ['% Generated by build_phase1_note.py; do not edit.']
    for key, value in v.items():
        lines.append(rf'\newcommand{{\V{key}}}{{{value}}}')
    (BUILD/'values.tex').write_text('\n'.join(lines)+'\n')

def compass(data, lang):
    t = TEXT[lang]
    fig, axes = plt.subplots(1, 3, figsize=(COLUMN, 1.45), subplot_kw=dict(projection='polar'))
    for ax, r, name in zip(axes, data['base'], t['modes']):
        for g in range(4):
            amp, ph = float(r[f'G{g+1}_amp']), np.deg2rad(float(r[f'G{g+1}_phase_deg']))
            ax.annotate('', xy=(ph, amp), xytext=(0, 0),
                        arrowprops=dict(arrowstyle='-|>', lw=.9, color=f'C{g}', shrinkA=0, shrinkB=0,
                                        mutation_scale=6))
        ax.set_ylim(0, 1.05); ax.set_yticks([.5, 1]); ax.set_yticklabels([])
        ax.set_xticks(np.deg2rad([0, 90, 180, 270])); ax.set_xticklabels(['0°', '', '180°', ''], fontsize=6)
        ax.tick_params(pad=-2); ax.grid(lw=.3)
        ax.set_title(f'{name}\n{float(r["frequency_Hz"]):.3f} Hz', pad=6)
    handles = [plt.Line2D([], [], color=f'C{g}', lw=1.2) for g in range(4)]
    fig.legend(handles, ['G1', 'G2', 'G3', 'G4'], ncol=4, loc='lower center', bbox_to_anchor=(.5, -.01),
               handlelength=1.2, columnspacing=1)
    fig.subplots_adjust(left=.07, right=.93, top=.74, bottom=.2, wspace=.62)
    fig.savefig(BUILD/f'fig_modes_{lang}.pdf'); plt.close(fig)

def fault_figure(data, lang):
    t, f, a = TEXT[lang], data['fault'], data['andes']
    fig, axes = plt.subplots(3, 1, figsize=(COLUMN, 2.45), sharex=True)
    for g in range(4):
        axes[0].plot(f['t_s'], f[f'G{g+1}_coi_deg'], label=f'G{g+1}')
    axes[0].legend(ncol=4, loc='upper right', handlelength=1.2, columnspacing=.9)
    axes[1].plot(f['t_s'], f['P7to8_MW'], color='k', label=t['ours'])
    axes[1].plot(a['t_s'], a['P7to8_MW'], color='C1', ls=(0, (3, 2)), lw=.9, label=t['andes'])
    axes[1].legend(loc='lower right', handlelength=1.8)
    coi = sum(H[g]*f[f'G{g+1}_df_Hz'] for g in range(4))/H.sum()
    coi_andes = sum(H[g]*a[f'G{g+1}_df_Hz'] for g in range(4))/H.sum()
    axes[2].plot(f['t_s'], coi, color='k')
    axes[2].plot(a['t_s'], coi_andes, color='C1', ls=(0, (3, 2)), lw=.9)
    for ax, label in zip(axes, [t['angle'], t['tie'], t['coi']]):
        ax.axvspan(*FAULT, color='0.85', lw=0); ax.set_ylabel(label, linespacing=1.1)
        ax.yaxis.set_label_coords(-.13, .5)
    axes[2].set_xlabel(t['time']); axes[2].set_xlim(0, 12)
    for k, ax in enumerate(axes):
        ax.text(.012, .95, f'({"abc"[k]})', transform=ax.transAxes, va='top', fontsize=7, zorder=5,
                bbox=dict(boxstyle='square,pad=.15', fc='white', ec='none'))
    fig.subplots_adjust(left=.17, right=.98, top=.99, bottom=.11, hspace=.12)
    fig.savefig(BUILD/f'fig_fault_{lang}.pdf'); plt.close(fig)

def sensitivity_figure(data, lang):
    t, rows_ = TEXT[lang], data['sens']
    fig, ax = plt.subplots(figsize=(COLUMN, 1.7))
    for k, (c, name) in enumerate(zip(data['cases'], t['cases'])):
        rr = [r for r in rows_ if r['sweep'] == 'transformer_x' and r['exciter'] == c['exciter']
              and int(r['tie_circuits']) == c['tie_circuits']]
        xs = [float(r['transformer_x_pu_machine']) for r in rr]
        ax.plot(xs, [float(r['inter_area_Hz']) for r in rr], color=f'C{k}', marker='o', ms=2.4,
                mfc='white', mew=.7, label=name)
        ax.axhline(c['ref2_Hz'], color=f'C{k}', lw=.6, ls=':')
        ax.plot([c['frequency_crossing_transformer_x_pu']], [c['ref2_Hz']], marker='D', ms=3.2,
                color=f'C{k}', ls='none')
        # Reference label right of the crossing, where the curve lies below the dotted line;
        # the lowest case is labelled under its line to stay clear of the next curve.
        below = k == 0
        ax.text(.104, c['ref2_Hz']+(-.003 if below else .003), f"{t['ref']}: {c['ref2_Hz']:.3f}",
                color=f'C{k}', fontsize=6.3, va='top' if below else 'bottom')
        ax.text(.153, float(rr[-1]['inter_area_Hz'])+(.012 if k == 1 else 0), name, color=f'C{k}', fontsize=6.3, va='center')
    ax.axvline(.15, color='0.4', lw=.6, ls='--')
    ax.text(.147, .395, t['textbook'], ha='right', va='bottom', fontsize=6.3, color='0.3')
    ax.set(xlabel=t['xt'], ylabel=t['freq'], xlim=(0, .215), ylim=(.39, .80))
    ax.set_xticks([0, .05, .1, .15])
    fig.subplots_adjust(left=.15, right=.98, top=.98, bottom=.19)
    fig.savefig(BUILD/f'fig_sensitivity_{lang}.pdf'); plt.close(fig)

def compile_note(lang):
    import pymupdf
    source = HERE/f'note_{lang}.tex'
    shutil.copy(source, BUILD/source.name)
    for _ in range(2):
        result = subprocess.run(['xelatex', '-interaction=nonstopmode', '-halt-on-error', source.name],
                                cwd=BUILD, capture_output=True, text=True)
        if result.returncode:
            raise RuntimeError(f'xelatex failed for {lang}; see {BUILD/source.with_suffix(".log").name}\n'
                               + result.stdout[-2000:])
    log = (BUILD/f'note_{lang}.log').read_text(errors='ignore')
    for warning in ('Overfull \\hbox', 'Reference `', 'Citation `', 'Missing character'):
        if warning in log:
            raise RuntimeError(f'[{lang}] LaTeX log reports: {warning}')
    target = PDF_OUT/f'Phase1_Note_{lang.upper()}.pdf'
    with pymupdf.open(BUILD/f'note_{lang}.pdf') as doc:
        if len(doc) != 2:  # check before publishing, so a failed build leaves no output
            raise RuntimeError(f'{target.name} has {len(doc)} pages; the plan asks for a two-page note.')
    shutil.copy(BUILD/f'note_{lang}.pdf', target)
    qa = ROOT/'qa'/target.stem
    qa.mkdir(parents=True, exist_ok=True)
    with pymupdf.open(target) as doc:
        for i, page in enumerate(doc, 1):
            page.get_pixmap(dpi=130).save(qa/f'{i}.png')
    print(target)

def main():
    BUILD.mkdir(exist_ok=True); PDF_OUT.mkdir(parents=True, exist_ok=True)
    data = collect()
    write_values(data['values'])
    for lang in ('en', 'zh'):
        compass(data, lang); fault_figure(data, lang); sensitivity_figure(data, lang)
        compile_note(lang)

if __name__ == '__main__':
    main()

"""Analyze actual runs; distinguish profiles and align event right limits explicitly."""
from pathlib import Path
import csv, json, hashlib, tarfile
import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'results'

def write(name, columns, rows):
    with (OUT/name).open('w',newline='') as stream:
        w=csv.writer(stream);w.writerow(columns);w.writerows(rows)

def main():
    historical=json.loads((OUT/'historical/observed.json').read_text())
    actual=sorted(historical['modes'],key=lambda m:m['freq'])
    current=np.genfromtxt(OUT/'official_193_modes.csv',delimiter=',',names=True)
    current=np.sort(current,order='frequency_Hz')
    published=[
        ('inter_area',-.19177,4.2247,.67239,4.5346,[5,4,5,4]),
        ('local_low',-.65252,6.8343,1.0877,9.5046,[5,4,4,4]),
        ('local_high',-.6564,7.086,1.1278,9.224,[4,3,4,3])]
    comparisons=[];rounded_checks=[]
    for ref,run,new in zip(published,actual,current):
        name,sigma,omega,freq,damping,precision=ref
        observed=[run[k] for k in ['sigma','omega','freq','damping']]
        passed=[abs(a-b)<=.5*10**(-p)+1e-12 for a,b,p in
                zip([sigma,omega,freq,damping],observed,precision)]
        rounded_checks.append(dict(mode=name,quantities=['sigma','omega','frequency','damping'],
                                   within_printed_rounding=passed))
        comparisons.append([name,sigma,run['sigma'],omega,run['omega'],freq,run['freq'],
            damping,run['damping'],new['frequency_Hz'],new['damping_pct'],
            (new['frequency_Hz']/run['freq']-1)*100,all(passed)])
    write('official_version_comparison.csv',
          ['mode','published_sigma','andes_175_sigma','published_omega','andes_175_omega',
           'published_frequency_Hz','andes_175_frequency_Hz','published_damping_pct',
           'andes_175_damping_pct','andes_193_frequency_Hz','andes_193_damping_pct',
           'version_frequency_change_pct','175_matches_printed_rounding'],comparisons)

    local_path=ROOT.parent/'CloudDRTS/results/detailed_fault.csv'
    local=np.genfromtxt(local_path,delimiter=',',names=True)
    groups=[('COI_angle','deg',[f'G{i}_coi_deg' for i in range(1,5)]),
            ('speed_deviation','Hz',[f'G{i}_df_Hz' for i in range(1,5)]),
            ('bus8_voltage','pu',['V8_pu']),('tie_power','MW',['P7to8_MW'])]
    metric_rows=[];side_rows=[]
    for tag in ['consistent','fine']:
        ext=np.genfromtxt(OUT/f'matched_fault_{tag}_waveform.csv',delimiter=',',names=True)
        record=json.loads((OUT/f'fault_{tag}_run.json').read_text())
        assert record['end_time']==12 and not record['busted']
        aligned=np.column_stack([np.interp(local['t_s'],ext['t_s'],ext[name])
                                 for name in local.dtype.names[1:]])
        # The external time-series saves pre-callback values at exactly tf/tc.
        # Replace only the two event rows by the independently solved right limits.
        # No sample is removed from the final error statistics.
        # Reconstruct the external initial absolute delta from its first event:
        # use recorded initial angles exported by run_fault.py.
        delta0=np.array(record['initial_delta_rad'])
        for event in record['events']:
            event_t=event['t'];index=np.where(np.isclose(local['t_s'],event_t,atol=1e-10,rtol=0))[0]
            assert len(index)==1
            right=event['right_limit']
            angles=np.rad2deg(np.array(right['delta_rad'])-delta0)
            angles-=np.average(angles,weights=[6.5,6.5,6.175,6.175])
            values=np.r_[angles,(np.array(right['omega_pu'])-1)*60,right['V8_pu'],right['P7to8_MW']]
            for j,name in enumerate(local.dtype.names[1:]):
                side_rows.append([tag,event_t,name,float(aligned[index[0],j]),float(values[j]),
                                  float(local[name][index[0]])])
            aligned[index[0]]=values
        errors=np.column_stack([aligned[:,i]-local[name] for i,name in enumerate(local.dtype.names[1:])])
        for label,unit,names in groups:
            ix=[list(local.dtype.names[1:]).index(n) for n in names]
            error=errors[:,ix]
            peak=float(np.max(abs(error)));rms=float(np.sqrt(np.mean(error**2)))
            metric_rows.append([tag,record['requested_step'],label,unit,len(local),peak,rms])
        write(f'aligned_fault_{tag}.csv',list(local.dtype.names),np.column_stack([local['t_s'],aligned]))
        write(f'fault_error_{tag}.csv',list(local.dtype.names),np.column_stack([local['t_s'],errors]))
    write('fault_waveform_comparison.csv',
          ['run','requested_step_s','quantity','unit','comparison_samples','max_abs_error','rms_error'],metric_rows)
    write('event_side_alignment.csv',['run','t_s','variable','raw_interpolated_left_value',
          'external_solved_right_value','local_right_value'],side_rows)
    xlsx_hash=hashlib.sha256((ROOT/'case/kundur_full.xlsx').read_bytes()).hexdigest()
    with tarfile.open(ROOT/'historical/andes-1.7.5.tar.gz') as archive:
        historical_hash=hashlib.sha256(archive.extractfile('andes-1.7.5/andes/cases/kundur/kundur_full.xlsx').read()).hexdigest()
    assert xlsx_hash==historical_hash
    checks=dict(
        historical_release=historical['andes'],
        published_output_version='1.7.5.post2.dev0+ga53a31fa (2022-07-05)',
        published_rounding_checks=rounded_checks,
        published_rounding_pass=all(all(c['within_printed_rounding']) for c in rounded_checks),
        official_xlsx_identical_across_versions=(xlsx_hash==historical_hash),
        official_xlsx_sha256=hashlib.sha256((ROOT/'case/kundur_full.xlsx').read_bytes()).hexdigest(),
        local_fault_csv_sha256=hashlib.sha256(local_path.read_bytes()).hexdigest(),
        same_parameter_modal=json.loads((OUT/'modal_comparison_checks.json').read_text()),
        waveform_alignment='Interpolate to all 2401 local samples; replace only exact t=1,1.08 rows by external right-limit event solves; omit no samples.',
        waveform_metrics=[dict(zip(['run','step_s','quantity','unit','samples','max_abs','rms'],r)) for r in metric_rows],
        default_andes_fault_failure=json.loads((OUT/'fault_run.json').read_text()),
        event_adapter='Cartesian current-balance event initialization, fixed differential states, checked against full ANDES g; fixed-field unsaturated GENROU active-I/reactive-Z profile only.',
        claim_limits=['Only experiment 1 phasor calculations receive this independent ANDES check.',
          'This is not the original 1991 paper configuration.',
          'The default ANDES deep-fault run failed; only the disclosed event-adapter runs completed.',
          'A factor-four time-domain error reduction is not established; version/model regularization remains.',
          'STEPSS published values are supplementary references and were not executed locally.'])
    (OUT/'verification_summary.json').write_text(json.dumps(checks,indent=2,ensure_ascii=False))
    assert checks['published_rounding_pass']
    print(json.dumps(dict(published_rounding_pass=checks['published_rounding_pass'],waveform_metrics=checks['waveform_metrics']),indent=2))

if __name__=='__main__':main()

"""Reproduce benchmarks with installed ANDES 1.9.3 and bundled 1.7.5 source."""
from pathlib import Path
import subprocess,sys,os,tarfile,hashlib,argparse
ROOT=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--rerun-local-transient',action='store_true',
                    help='also regenerate the custom DOP853 waveform; otherwise use supplied v5 data')
args=parser.parse_args()
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',MPLCONFIGDIR=str(ROOT/'matplotlib'))
(ROOT/'logs').mkdir(exist_ok=True)

def run(script,arguments=(),extra_env=None):
    log=ROOT/'logs'/f'{Path(script).stem}_{"_".join(arguments).replace("/","_") or "run"}.log'
    with log.open('w') as stream:
        subprocess.run([sys.executable,'-u',str(ROOT/script),*arguments],env=extra_env or env,
                       stdout=stream,stderr=subprocess.STDOUT,check=True)
    print(f'Completed {Path(script).name}: {log.name}',flush=True)

if args.rerun_local_transient:
    run('../CloudDRTS/detailed_transient.py')
archive=ROOT/'historical/andes-1.7.5.tar.gz'
assert hashlib.sha256(archive.read_bytes()).hexdigest()=='878e187db89484f6e778ea19f31a8830440f7bef28d4a59e115847d54092b678'
source=ROOT/'historical/andes-1.7.5'
if not source.exists():
    with tarfile.open(archive) as bundle:bundle.extractall(ROOT/'historical',filter='data')
run('run_historical.py',extra_env=dict(env,PYTHONPATH=str(source)))
run('run_comparison.py')
run('run_fault.py')
run('run_fault.py',['--step','0.0005','--tag','fine'])
# analyze_results.py records the default (unadapted) ANDES fault run, which is
# expected to stop at the stability criterion. Run it here so a clean clone has
# results/fault_run.json; a non-zero exit is the documented outcome.
failed_log=ROOT/'logs'/'run_fault_default_failed_run.log'
with failed_log.open('w') as stream:
    code=subprocess.run([sys.executable,'-u',str(ROOT/'run_fault_default_failed.py')],env=env,
                        stdout=stream,stderr=subprocess.STDOUT).returncode
if not (ROOT/'results'/'fault_run.json').exists():
    raise RuntimeError(f'Default fault record missing; see {failed_log}')
print(f'Completed run_fault_default_failed.py (exit {code}, failure expected): {failed_log.name}',flush=True)
run('analyze_results.py')
run('plot_results.py')
print('All requested benchmark calculations and numerical analysis completed.')

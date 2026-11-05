"""Run one native VPSC-DDRX calculation per rate, then plot the saved states."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
import json
import os
from pathlib import Path
import subprocess
import shutil
import sys
import time
import tempfile

os.environ.setdefault('OMP_NUM_THREADS','1')
os.environ['OPENBLAS_NUM_THREADS']='1'
from app.model import ROOT, require, save_json, sha, validate, write_case
from app.outputs import FIGURES, DATA, export_outputs, output_choices
from app.stock import uses_original
import numpy as np

def table(path):
    a=np.genfromtxt(path,delimiter=',',names=True,deletechars='',encoding='ascii',ndmin=1)
    require(len(a)>0 and all(np.isfinite(a[k]).all() for k in a.dtype.names),f'Incomplete or nonfinite output in {path}.')
    return a

def check(c,case,steps,rate):
    history=table(case/'MG_HISTORY.csv')
    mode=table(case/'MG_MODE_HISTORY.csv')
    act=table(case/'SYSTEM_ACTIVITY.csv')
    drx=table(case/'DRX_HISTORY.csv')
    e=c['experiment']
    strain=history['axial_strain']
    for series in (history,mode,act,drx):
        require(len(series)==steps+1,'Solver stopped before the requested endpoint. Read console.log.')
    require(np.all(np.diff(strain)>0) and abs(strain[-1]-e['end_strain'])<1e-8,'Incorrect axial strain endpoint.')
    require(np.allclose(mode['axial_strain'],strain) and np.allclose(act['axial_strain'],strain) and np.allclose(drx['time_s']*rate,strain),'Native histories use different time points.')
    require(min(drx['drx_fraction'])>=0 and max(drx['drx_fraction'])<=1 and (c['drx']['enabled'] or max(drx['drx_fraction'])<1e-12),'Invalid DRX volume fraction.')
    first=1
    for i,m in enumerate(c['slip_modes'],1):
        last=first+len(m['systems'])
        shares=sum(act[f'system_{j}'] for j in range(first,last))
        require(np.allclose(shares,history[f'family_{i}'],atol=1e-10,rtol=0),f'Slip-system activity does not sum to family {i}.')
        first=last
    require(np.allclose(sum(act[k] for k in act.dtype.names[1:]),1,atol=1e-10,rtol=0),'Slip activity does not sum to one.')
    require(np.allclose(sum(mode[f'rho_{i}'] for i in range(1,len(c['slip_modes'])+1)),drx['rho_mean_m-2'],rtol=1e-10),'Mode densities do not sum to total density.')
    import csv
    with (case/'MG_SNAPSHOTS.csv').open() as f:
        maps=list(csv.DictReader(f))
    targets=c['plots']['strains']
    require(len(maps)==len(targets),'Missing requested strain snapshots.')
    for t,row in zip(targets,maps):
        actual=float(row['actual_axial_strain'])
        require(abs(float(row['target_strain'])-t)<1e-9 and -1e-9<=actual-t<=e['end_strain']/steps+1e-8 and (case/row['filename'].strip()).is_file(),'Snapshot does not match requested strain.')
    return dict(steps=steps,rate_per_s=rate,final_stress_MPa=float(history['axial_stress_MPa'][-1]),
                peak_stress_MPa=float(history['axial_stress_MPa'].max()),final_DRX_fraction=float(drx['drx_fraction'][-1]),
                final_mean_diameter_um=float(drx['mean_diameter_volume_um'][-1]),final_density_per_m2=float(drx['rho_mean_m-2'][-1]))

def run_case(c,rate,out,stock_path=None):
    case=out/f'rate_{rate}'
    case.mkdir()
    steps=write_case(c,rate,case,stock_path)
    print(f'Running {rate:g}/s for {steps} increments.',flush=True)
    start=time.monotonic()
    if not uses_original(c):
        with (case/'console.log').open('w') as log:
            proc=subprocess.run([str(ROOT/'bin/vpsc_drx.exe')],cwd=case,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL)
        require(proc.returncode==0 and (case/'MG_SNAPSHOTS.csv').is_file(),f'Native solver failed. See {case/"console.log"}')
        result=check(c,case,steps,rate)
    else:
        from app.stock import execute_stock
        result=execute_stock(c,case,steps,rate)
    result['elapsed_seconds']=time.monotonic()-start
    save_json(case/'result.json',result)
    print(f'Finished {rate:g}/s. Final axial stress {result["final_stress_MPa"]:.3f} MPa.',flush=True)
    return result

def completed_settings(out):
    c=json.loads((out/'settings_used.json').read_text(encoding='utf8'))
    first_rate=c['experiment']['rates_per_s'][0]
    meta=json.loads((out/f'rate_{first_rate}'/'effective_inputs.json').read_text(encoding='utf8'))
    if not c['drx']['enabled'] and meta.get('engine')!='original_vpsc':
        # Saved runs from the former population-only switch retain their own outputs.
        validate(dict(c,drx=dict(c['drx'],enabled=True)))
    else:
        validate(c,meta.get('material'))
    for rate in c['experiment']['rates_per_s']:
        require((out/f'rate_{rate}'/'result.json').is_file(),f'Missing completed run at {rate}/s.')
    return c


def plot_run(out):
    c = completed_settings(out)
    from app.plots import render
    render(c, out, out, set(output_choices(c)['figures']), set())


def export_saved(run, dest, options):
    c = completed_settings(run)
    c['outputs'] = options
    output_choices(c)
    require(not dest.exists(), f'Choose a new export folder: {dest}')
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Publish a complete export, including its manifest, in one rename.
    with tempfile.TemporaryDirectory(dir=dest.parent) as tmp:
        staging = Path(tmp) / 'export'
        staging.mkdir()
        export_outputs(c, run, staging)
        staging.rename(dest)
    print(f'Exported to {dest}', flush=True)

def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--settings',type=Path,default=ROOT/'settings.json')
    a.add_argument('--output',type=Path,help='New result folder; existing folders are preserved.')
    post = a.add_mutually_exclusive_group()
    post.add_argument('--plots-only',type=Path,metavar='COMPLETED_RUN')
    post.add_argument('--export-from',type=Path,metavar='COMPLETED_RUN', help='Export saved solver records without rerunning.')
    a.add_argument('--figures', nargs='*', choices=list(FIGURES), help='Figure selection for --export-from (empty selects none).')
    a.add_argument('--data', nargs='*', choices=list(DATA), help='Numerical output selection for --export-from.')
    a.add_argument('--workers',type=int,default=1)
    a.add_argument('--validate',action='store_true')
    args=a.parse_args()
    if args.export_from:
        run = args.export_from.resolve()
        options = dict(output_choices(json.loads((run/'settings_used.json').read_text(encoding='utf8'))))
        for key in ('figures','data'):
            if getattr(args, key) is not None:
                options[key] = getattr(args, key)
        dest = (args.output or run/'exports'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')).resolve()
        export_saved(run, dest, options)
        return
    require(args.figures is None and args.data is None, 'Use --figures and --data with --export-from, or edit outputs in the run settings.')
    if args.plots_only:
        plot_run(args.plots_only.resolve())
        return
    c=json.loads(args.settings.read_text(encoding='utf-8-sig'))
    validate(c)
    from app.model import grains
    grains(c['experiment'])
    if args.validate:
        print('Settings and input grains passed validation.')
        return
    require(args.workers>0,'Workers must be positive.')
    out=(args.output or ROOT/'results'/datetime.now().strftime('%Y%m%d_%H%M%S_%f')).resolve()
    out.mkdir(parents=True,exist_ok=False)
    save_json(out/'settings_used.json',c)
    print(out,flush=True)
    engine='original_vpsc' if uses_original(c) else 'vpsc_drx'
    binary='vpsc_original.exe' if uses_original(c) else 'vpsc_drx.exe'
    stock_path=None
    if uses_original(c):
        from app.stock import material_path
        stock_path=out/'original_material.sx'
        shutil.copyfile(material_path(c),stock_path)
    log=dict(status='running',engine=engine,native_executable_sha256=sha(ROOT/'bin'/binary),settings_sha256=sha(out/'settings_used.json'))
    save_json(out/'run_status.json',log)
    try:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            log['cases']=list(pool.map(lambda rate: run_case(c,rate,out,stock_path),c['experiment']['rates_per_s']))
        log['status']='solver_complete'
        save_json(out/'run_status.json',log)
        export_outputs(completed_settings(out), out, out)
        log['status']='complete'
    except BaseException as exc:
        log.update(status='failed',error=str(exc))
        raise
    finally:
        save_json(out/'run_status.json',log)
    print(f'Outputs: {out}',flush=True)

if __name__=='__main__':
    try:
        main()
    except Exception as exc:
        print(f'ERROR: {exc}',file=sys.stderr,flush=True)
        sys.exit(1)


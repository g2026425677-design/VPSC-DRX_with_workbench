"""Stock VPSC uses its own byte-for-byte material deck and its native hardening law."""
import shutil
import re
import subprocess
import tempfile
from pathlib import Path
import csv
import numpy as np
from .model import ROOT, require, sha, save_json

DEFAULT_SX = 'data/original_vpsc/Mg.sx'
ORIGINAL = ROOT/'bin/vpsc_original.exe'
REPORTER = ROOT/'bin/vpsc_original_report.exe'

def uses_original(c):
    return not c['drx']['enabled'] and c['drx'].get('hardening_when_off','original_material')!='paper_km'

def material_path(c):
    value=c['drx'].get('original_material_file',DEFAULT_SX)
    p=(ROOT/Path(value)).resolve()
    require(p.is_file(),f'Original VPSC material file does not exist: {p}')
    return p

def selected_material(c):
    info=inspect_material(material_path(c))
    expected={'original_dd_rev':23,'original_voce':0}.get(c['drx'].get('hardening_when_off'))
    require(expected is None or info['hardening_law']==expected,'The selected hardening law and .sx file disagree. Select a matching file or choose original_material.')
    return info

def native_table(path):
    def number(value):
        # Legacy Fortran E fields may omit E for three-digit exponents.
        return float(re.sub(r'^([+-]?\d*\.\d+)([+-]\d{3})$',r'\1E\2',value.strip()))
    return np.loadtxt(path,skiprows=1,ndmin=2,encoding='ascii',converters=number)

def inspect_material(path):
    # The Fortran reader remains authoritative. This reads the same header to
    # label native plots and check the one-phase workbench and plot conventions.
    lines=iter(path.read_text(encoding='latin1').splitlines())
    def nums(n,cast=float):
        return [cast(s.replace('D','E').replace('d','e')) for s in next(lines).split()[:n]]
    name=next(lines).strip()
    token=next(lines).strip().lower()[:5]
    system={'cubic':'cubic','hexag':'hexagonal','trigo':'trigonal','tetra':'tetragonal','ortho':'orthorhombic'}.get(token)
    require(system is not None,'Select an original cubic, hexagonal, trigonal, tetragonal or orthorhombic material.')
    cell=nums(6)
    require(len(cell)==6,'Invalid original unit cell.')
    nind=4 if system in ('hexagonal','trigonal') else 3
    for _ in range(10): next(lines)  # stiffness label + 6 rows, expansion label + row, modes label
    total=nums(1,int)[0]
    active_count=nums(1,int)[0]
    active=nums(active_count,int)
    modes=[]
    for i in range(1,total+1):
        label=next(lines).strip()
        mode,count,sense,twin=nums(4,int)
        require(mode==i,'Original VPSC mode IDs must be consecutive.')
        if twin: next(lines)
        systems=[nums(2*nind,int) for _ in range(count)]
        if mode in active:
            modes.append(dict(name=label,systems=systems,twin=bool(twin)))
    require(len(modes)==active_count,'Missing active mode in the original material.')
    next(lines)
    law=nums(1,int)[0]
    nums(1,int)  # IRATESENS
    require(law in (0,1,20,23),'The workbench supports native single-phase Voce, MTS, DD and DD_REV material files.')
    diameter=None
    if law==23:
        next(lines)  # Backstress factors
        next(lines)  # qBS, mREV
        diameter=nums(1)[0]*1e6
    return dict(name=name,system=system,lattice_lengths_angstrom=cell[:3],lattice_angles_deg=cell[3:],
                modes=modes,hardening_law=law,hardening_grain_size_um=diameter)

def execute_stock(c,case,steps,rate):
    with (case/'console.log').open('w') as f:
        proc=subprocess.run([str(ORIGINAL)],cwd=case,stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL)
    require(proc.returncode==0 and (case/'STR_STR.OUT').is_file(),f'Original VPSC failed. See {case/"console.log"}')
    return collect_stock(c,case,steps,rate)

def collect_stock(c,case,steps,rate):
    import json
    meta=json.loads((case/'effective_inputs.json').read_text(encoding='utf8'))
    material=meta['material']
    native=native_table(case/'STR_STR.OUT')
    require(len(native)==steps+1 and np.isfinite(native).all() and abs(-native[-1,4]-c['experiment']['end_strain'])<=max(1e-8,c['experiment']['end_strain']*6e-5),
            f'Original VPSC did not reach the requested endpoint. See {case/"console.log"}')
    # Keep all authoritative stock files. The replay exports additional full-
    # precision observations only if its standard stress, activity and texture agree.
    checks={}
    with tempfile.TemporaryDirectory(dir=case) as td:
        report=Path(td)
        for name in ('VPSC8.in','compression.dat','initial.tex','magnesium.sx','stock_snapshots.in'):
            shutil.copyfile(case/name,report/name)
        with (case/'reporter.log').open('w') as f:
            proc=subprocess.run([str(REPORTER)],cwd=report,stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL)
        require(proc.returncode==0 and (report/'MG_SNAPSHOTS.csv').is_file(),f'Original reporting replay failed. See {case/"reporter.log"}')
        names=['STR_STR.OUT','ACT_PH1.OUT','TEX_PH1.OUT']
        if (case/'DIS_DEN.OUT').is_file(): names.append('DIS_DEN.OUT')
        for name in names:
            checks[name]=sha(case/name)
            require(checks[name]==sha(report/name),f'Reporting replay differs from the original {name}. Original output retained.')
        for name in ['MG_HISTORY.csv','SYSTEM_ACTIVITY.csv','MG_SNAPSHOTS.csv']:
            shutil.copyfile(report/name,case/name)
        for p in [*report.glob('MG_STATE_*.csv'),*report.glob('STOCK_DENSITY.csv')]:
            shutil.copyfile(p,case/p.name)
    hist=np.loadtxt(case/'MG_HISTORY.csv',delimiter=',',skiprows=1,ndmin=2)
    activity=np.loadtxt(case/'SYSTEM_ACTIVITY.csv',delimiter=',',skiprows=1,ndmin=2)
    # STR_STR.OUT displays 5 significant digits. Full-precision reporter records
    # come from the same solve and must agree to that native display precision.
    require(len(hist)==steps+1 and np.isfinite(hist).all() and np.allclose(hist[:,1],-native[:,10],rtol=5.1e-5,atol=1e-10),'Axial stress differs from stock STR_STR.OUT.')
    require(len(activity)==len(hist) and np.allclose(activity[:,0],hist[:,0]) and np.allclose(activity[:,1:].sum(axis=1),1),'Invalid original system activity.')
    start=1
    for i,mode in enumerate(material['modes'],2):
        end=start+len(mode['systems'])
        require(np.allclose(hist[:,i],activity[:,start:end].sum(axis=1),atol=1e-10,rtol=0),'System and native family activities disagree.')
        start=end
    with (case/'MG_SNAPSHOTS.csv').open() as f:
        snapshots=list(csv.DictReader(f))
    targets=c['plots']['strains']
    require(len(targets)==len(snapshots),'Missing stock orientation snapshots.')
    for row,target in zip(snapshots,targets):
        actual=float(row['actual_axial_strain'])
        require(abs(float(row['target_strain'])-target)<1e-9 and -1e-9<=actual-target<=meta['increment']+1e-8
                and (case/row['filename'].strip()).is_file(),'Stock orientation snapshot missed its strain.')
    rho=None
    if material['hardening_law'] in (20,23):
        density=np.loadtxt(case/'STOCK_DENSITY.csv',delimiter=',',skiprows=1,ndmin=2)
        require(len(density)==len(hist) and np.isfinite(density).all() and np.allclose(density[:,0],hist[:,0]),'Invalid stock density history.')
        rho=float(density[-1,-1])
        if material['hardening_law']==23:
            ref=native_table(case/'DIS_DEN.OUT')
            # Native file contains post-update states only and uses E14.4.
            require(len(ref)==steps and np.allclose(density[1:,-1],ref[:,5],rtol=5.2e-4,atol=0),'DD_REV density disagrees with native DIS_DEN.OUT.')
    save_json(case/'stock_parity.json',dict(stock_sha256=sha(ORIGINAL),reporter_sha256=sha(REPORTER),
              material_sha256=sha(case/'magnesium.sx'),matched_native_outputs=checks,records=len(native)))
    return dict(engine='original_vpsc',steps=steps,rate_per_s=rate,peak_stress_MPa=float(hist[:,1].max()),
                final_stress_MPa=float(hist[-1,1]),final_DRX_fraction=0.0,
                final_mean_diameter_um=None,hardening_grain_size_um=material['hardening_grain_size_um'],final_density_per_m2=rho)


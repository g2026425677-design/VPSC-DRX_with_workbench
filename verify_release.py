"""Recompare the delivered example with the frozen calibration and exercise editable controls."""
import argparse
import copy
import io
import json
import tempfile
import zipfile
import numpy as np
from simulate import ROOT, run_case, plot_run
from app.model import save_json, sha, validate, require

def compare(out):
    baseline=ROOT/'pre_hollow_calibration/original_delivery.zip'
    fingerprint='427fafe32cd8c8e8acf39476a2e74a2984816796e0c1908f8ce01f44ecd90536'
    require(sha(baseline)==fingerprint,'Frozen delivery has changed.')
    c=json.loads((out/'settings_used.json').read_text(encoding='utf8'))
    default=json.loads((ROOT/'app/default_settings.json').read_text(encoding='utf8'))
    for obj in (c,default):
        obj.pop('outputs',None)  # Presentation settings do not enter the archive comparison.
        obj['drx'].pop('hardening_when_off',None)
        obj['drx'].pop('original_material_file',None)
    require(c==default,'Compare the original calibrated settings only.')
    rows=[]
    with zipfile.ZipFile(baseline) as z:
        for rate in c['experiment']['rates_per_s']:
            case=out/f'rate_{rate}'
            pre='VPSC_code-main/mg_calibration/paper_revision_20261103/validation/auxiliary_final/'+case.name+'/'
            errors={}
            for name in ['MG_HISTORY.csv','DRX_HISTORY.csv','MG_MODE_HISTORY.csv','DRX_CYCLES.csv',
                         *[f'MG_STATE_{i:03d}.csv' for i in range(1,6)]]:
                old=np.loadtxt(io.BytesIO(z.read(pre+name)),delimiter=',',skiprows=1)
                new=np.loadtxt(case/name,delimiter=',',skiprows=1)
                require(old.shape==new.shape,f'{case.name}/{name} shape mismatch.')
                errors[name]=float(np.abs(new-old).max())
                require(np.array_equal(old,new),f'{case.name}/{name} differs from archive.')
            rows.append(dict(rate_per_s=rate,maximum_absolute_errors=errors,history_records=len(np.loadtxt(case/'MG_HISTORY.csv',delimiter=',',skiprows=1))))
    return rows

def controls():
    from app.plots import mt
    p=json.loads((ROOT/'app/default_settings.json').read_text(encoding='utf8'))
    results={}
    with tempfile.TemporaryDirectory(dir=ROOT/'verification') as td:
        from pathlib import Path
        for name in ['archived_Mg','changed_Mg','cubic_FCC']:
            case=Path(td)/name
            case.mkdir()
            c=copy.deepcopy(p)
            c['experiment'].update(end_strain=.003,rates_per_s=[.2])
            c['plots'].update(strains=[0,.00135,.003])
            if name!='archived_Mg':
                c['drx']['enabled']=False
                c['drx']['hardening_when_off']='paper_km'
                c['slip_modes'][0]['friction_MPa']*=1.3
                c['experiment']['orientation_convention']='aztec'
                c['plots']['pole_planes'].append([1,0,-1,0])
            if name=='cubic_FCC':
                c['crystal'].update(system='cubic',lattice_lengths_angstrom=[3.6]*3,lattice_angles_deg=[90]*3)
                c['plots']['pole_planes']=[[1,1,1]]
                c['experiment']['initial_grains']='equal_spheres_number_mean'
                m=c['slip_modes'][0]
                m['systems']=[[1,1,1,0,1,-1],[1,1,1,1,0,-1],[1,1,1,1,-1,0],[-1,1,1,0,1,-1],
                     [-1,1,1,1,0,1],[-1,1,1,1,1,0],[-1,-1,1,0,1,1],[-1,-1,1,1,0,1],
                     [-1,-1,1,1,-1,0],[1,-1,1,0,1,1],[1,-1,1,1,0,-1],[1,-1,1,1,1,0]]
                m['name']='FCC {111}<110>'
                c['slip_modes']=[m]
            validate(c)
            save_json(case/'settings_used.json',c)
            record=run_case(c,.2,case)
            raw=json.loads((case/'rate_0.2/effective_inputs.json').read_text())['native_namelist']
            assert raw['enabled']==int(c['drx']['enabled']) and raw['mode_tau'][0]==c['slip_modes'][0]['friction_MPa']
            if name!='archived_Mg':
                plot_run(case)
                assert (case/'figures/texture_rate_0.2_state_03.png').is_file()
                record.update(mtex_version=mt.__version__,generated_figures=len(list((case/'figures').glob('*.png'))))
            if name=='cubic_FCC':
                grains=np.loadtxt(case/'rate_0.2/MG_INITIAL_GRAINS.csv',delimiter=',',skiprows=1)
                assert np.ptp(grains[:,1])==0 and np.allclose(grains[:,2],1/328) and np.allclose(grains[:,3],1)
            results[name]=record
        require(abs(results['archived_Mg']['final_stress_MPa']-results['changed_Mg']['final_stress_MPa'])>.01,'Edited Mg test did not change stress.')
    return results

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--controls-only',action='store_true')
    a=ap.parse_args()
    dest=ROOT/'verification'
    dest.mkdir(exist_ok=True)
    save_json(dest/'editable_controls.json',controls())
    if not a.controls_only:
        save_json(dest/'archived_comparison.json',compare(ROOT/'results/calibrated_example'))
    print('Verification passed.')


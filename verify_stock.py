"""Verify the original-VPSC DRX-off branch and the user-selectable hardening laws."""
from copy import deepcopy
import json
from pathlib import Path
import tempfile
import numpy as np
from app.model import ROOT, validate, sha, save_json, write_case
from app.stock import uses_original
from simulate import run_case, plot_run

def full_comparison():
    out=ROOT/'results/original_vpsc_off'
    direct=ROOT/'verification/original_vpsc/direct'
    c=json.loads((out/'settings_used.json').read_text())
    original=ROOT.parent/'VPSC_code-main (1)'/'VPSC_code-main (1)'/'VPSC_code-main'
    files={}
    for name in ('vpsc8dim.for','vpsc8main.for','vpsc8sub.for'):
        a=sha(original/name)
        assert a==sha(ROOT/'source/stock_vpsc'/name)
        files[name]=a
    sx=original/'ex13_dd_rev/Mg/Mg.sx'
    records=[]
    for rate in c['experiment']['rates_per_s']:
        run=out/f'rate_{rate}'
        ref=direct/f'rate_{rate}'
        assert sha(run/'magnesium.sx')==sha(sx)
        matches={}
        for name in ['initial.tex','compression.dat','VPSC8.in','magnesium.sx',
                     'STR_STR.OUT','ACT_PH1.OUT','TEX_PH1.OUT','DIS_DEN.OUT']:
            matches[name]=sha(run/name)
            assert matches[name]==sha(ref/name),(rate,name)
        records.append(dict(rate_per_s=rate,matched_hashes=matches,
                            result=json.loads((run/'result.json').read_text())))
    return dict(source_matches_requested_original_folder=files,material_sha256=sha(sx),cases=records)

def control_checks():
    base=json.loads((ROOT/'settings.json').read_text())
    records={}
    with tempfile.TemporaryDirectory(dir=ROOT/'verification') as tmp:
        tmp=Path(tmp)
        text=(ROOT/'data/original_vpsc/Mg.sx').read_bytes().decode('latin1')
        edited=tmp/'edited.sx'
        edited.write_bytes(text.replace('118.0','160.0',1).encode('latin1'))
        for label in ('original_dd_rev','edited_dd_rev','original_voce','paper_km'):
            c=deepcopy(base)
            c['drx'].update(enabled=False,hardening_when_off='original_dd_rev')
            c['experiment'].update(rates_per_s=[.2],end_strain=.003)
            c['plots']['strains']=[0,.00135,.003]
            if label=='edited_dd_rev': c['drx']['original_material_file']=str(edited)
            if label=='original_voce':
                c['drx'].update(hardening_when_off=label,original_material_file='data/original_vpsc/Mg_voce.sx')
            if label=='paper_km': c['drx']['hardening_when_off']=label
            validate(c)
            case=tmp/label
            case.mkdir()
            save_json(case/'settings_used.json',c)
            result=run_case(c,.2,case)
            plot_run(case)
            assert (case/'figures/texture_rate_0.2_state_03.png').is_file()
            run=case/'rate_0.2'
            assert (run/'stock_parity.json').exists()==uses_original(c)
            assert (run/'drx.nml').exists()==(label=='paper_km')
            result['figure_count']=len(list((case/'figures').glob('*.png')))
            records[label]=result
        assert abs(records['edited_dd_rev']['final_stress_MPa']-records['original_dd_rev']['final_stress_MPa'])>1
        # Turning DRX on still generates byte-identical physics inputs for the calibrated executable.
        c=deepcopy(base)
        expected=ROOT/'results/calibrated_example'
        for rate in c['experiment']['rates_per_s']:
            case=tmp/f'on_{rate}'
            case.mkdir()
            write_case(c,rate,case)
            for name in ('VPSC8.in','compression.dat','initial.tex','initial_sizes.dat','drx.nml','magnesium.sx'):
                assert (case/name).read_bytes()==(expected/f'rate_{rate}'/name).read_bytes(),(rate,name)
    return records

if __name__=='__main__':
    report=dict(full_original_comparison=full_comparison(),hardening_control_checks=control_checks())
    save_json(ROOT/'verification/original_vpsc/verification.json',report)
    print('Original-source parity and user-defined hardening checks passed.')


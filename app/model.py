"""Translate one editable parameter file into the native single-phase compression case."""
import csv
import hashlib
import json
import math
import shutil
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TEMPLATES = Path(__file__).resolve().parent / 'templates'
SYMMETRY = {'hexagonal':'6/mmm', 'trigonal':'-3m', 'cubic':'m-3m',
            'tetragonal':'4/mmm', 'orthorhombic':'mmm'}

def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False)+'\n', encoding='utf-8')

def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()

def require(test, message):
    if not test:
        raise ValueError(message)

def indexed(v, n):
    return isinstance(v,list) and len(v)==n and all(type(i) is int for i in v) and any(v) and (n!=4 or sum(v[:3])==0)

def canonical(v):
    sign=1 if next(i for i in v if i)>0 else -1
    divisor=math.gcd(*v)
    return tuple(sign*i//divisor for i in v)

def validate(c,stock_material=None):
    schema=json.loads((ROOT/'app/default_settings.json').read_text(encoding='utf8'))
    from .outputs import output_choices
    output_choices(c)
    require(c.keys()-{'outputs'}==schema.keys()-{'outputs'}, 'Unknown or missing section in settings.')
    for k in c.keys()-{'slip_modes','outputs'}:
        optional={'original_material_file','hardening_when_off'} if k=='drx' else set()
        require(c[k].keys()-optional==schema[k].keys()-optional, f'Unknown or missing field in {k}.')
    x,e,d,p=c['crystal'],c['experiment'],c['drx'],c['plots']
    require(type(d['enabled']) is bool,'DRX enabled must be true or false.')
    from .stock import uses_original, selected_material
    require(d.get('hardening_when_off','original_material') in ('original_material','original_dd_rev','original_voce','paper_km'),'Unknown hardening law for DRX-off calculations.')
    if not uses_original(c):
        _validate_km(c)
        system=x['system']
    else:
        stock_material=stock_material or selected_material(c)
        system=stock_material['system']
    n=4 if system in ('hexagonal','trigonal') else 3
    require(e['temperature_C']>-273.15,'Temperature must exceed absolute zero.')
    rates=e['rates_per_s']
    require(len(rates)>0 and min(rates)>0 and len(set(rates))==len(rates),'Use unique, positive axial strain rates.')
    require(e['end_strain']>0 and 0<e['max_strain_increment']<=e['end_strain'] and type(e['seed']) is int and 1<=e['seed']<=2147483646,'Invalid strain, increment or seed.')
    require(e['orientation_convention'] in ('archived','aztec','vpsc') and e['initial_grains'] in ('measured_cohorts','equal_spheres_number_mean','equal_spheres_area_mean'),'Unknown orientation or grain representation.')
    ss=p['strains']
    require(1<=len(ss)<=16 and ss==sorted(set(ss)) and min(ss)>=0 and max(ss)<=e['end_strain'],'Use 1–16 ordered, distinct snapshot strains in the simulated interval.')
    require(0<p['halfwidth_deg']<=90 and p['pole_planes'] and all(indexed(v,n) for v in p['pole_planes']),'Invalid texture smoothing or pole plane.')
    require(p['ipf_specimen_directions'] and all(len(v)==3 and np.isfinite(v).all() and np.linalg.norm(v)>0 for v in p['ipf_specimen_directions']),'IPFs require nonzero specimen directions.')
    # Reject nan/inf in every numeric field, including fields not compared above.
    json.dumps(c,allow_nan=False)

def _validate_km(c):
    x,d,e=c['crystal'],c['drx'],c['experiment']
    schema=json.loads((ROOT/'app/default_settings.json').read_text(encoding='utf8'))
    system=x['system']
    require(isinstance(x['name'],str) and x['name'].isascii() and x['name'].strip() and '\n' not in x['name'],'Use a single-line ASCII material name in the native crystal file.')
    require(system in SYMMETRY, f'Available crystal types: {list(SYMMETRY)}')
    hcp=system in ('hexagonal','trigonal')
    n=4 if hcp else 3
    lengths=np.asarray(x['lattice_lengths_angstrom'],dtype=float)
    angles=x['lattice_angles_deg']
    require(lengths.shape==(3,) and np.isfinite(lengths).all() and min(lengths)>0,'Supply three positive lattice lengths in Å.')
    require(np.array_equal(angles,[90,90,120 if hcp else 90]),'Lattice angles do not match the selected crystal setting.')
    if hcp or system in ('cubic','tetragonal'):
        require(np.isclose(lengths[0],lengths[1],rtol=1e-9),'This crystal setting requires a=b.')
    if system=='cubic':
        require(np.isclose(lengths[0],lengths[2],rtol=1e-9),'Cubic symmetry requires a=b=c.')
    require(-1<x['poisson_ratio']<.5 and min(x['mu_reference_Pa'],x['chi'],x['burgers_reference_m'],x['recovery_reference_rate_per_s'])>0
            and type(x['slip_exponent']) is int and x['slip_exponent']>=2,'Invalid elastic, forest or slip-law parameter.')
    modes=c['slip_modes']
    require(1<=len(modes)<=8 and sum(len(m['systems']) for m in modes)<=48,'Native capacity is 8 modes and 48 systems.')
    seen=set()
    for m in modes:
        require(m.keys()==schema['slip_modes'][0].keys(),'Unknown or missing slip-mode field.')
        require(isinstance(m['name'],str) and m['name'].isascii() and m['name'].strip() and '\n' not in m['name'] and bool(m['systems']),'Supply an ASCII mode label and at least one slip system.')
        require(m['friction_MPa']>=0 and m['recovery_k2']>=0 and min(m['chi_ratio'],m['storage_k1_per_m'],m['initial_density_per_m2'],m['nucleus_density_per_m2'])>0,'Invalid mode hardening or density.')
        b=m['burgers_ratio']
        require((b=='hcp_c_plus_a' and hcp) or (isinstance(b,(float,int)) and b>0),'Burgers ratio must be positive or hcp_c_plus_a for a hexagonal cell.')
        for row in m['systems']:
            plane,slip=row[:n],row[n:]
            require(indexed(plane,n) and indexed(slip,n) and np.dot(plane,slip)==0,f'Invalid plane/direction or zone law in {m["name"]}: {row}')
            key=(canonical(plane),canonical(slip))
            require(key not in seen,f'Duplicate bidirectional slip system: {row}')
            seen.add(key)
    require(d['orientation_rule'] in ('gaussian','random') and 0<=d['orientation_spread_deg']<=180,'Invalid DRX orientation rule.')
    require(min(d['threshold_density_per_m2'],d['nucleus_diameter_m'],d['line_energy_factor'])>0 and min(d['mobility_m4_per_J_s'],d['nucleation_coefficient'])>=0,'Invalid DRX coefficient.')
    require(e['temperature_C']>-273.15 and shear_modulus(c)>0,'Nonphysical temperature or shear modulus.')

def shear_modulus(c):
    x=c['crystal']
    return x['mu_reference_Pa']*(1-x['mu_slope_per_K']*(c['experiment']['temperature_C']+273.15-x['mu_reference_temperature_K']))

def grains(e):
    path=ROOT/Path(e['grain_csv'])
    with path.open(encoding='utf-8-sig',newline='') as f:
        rows=list(csv.DictReader(f))
    if rows and rows[0]['Mean Orientation Euler1'].strip().lower() in ('°','deg','degrees',''):
        rows=rows[1:]
    angles=np.array([[float(r[f'Mean Orientation Euler{i}']) for i in (1,2,3)] for r in rows])
    sizes=np.array([float(r['Equivalent Circle Diameter'])*1e-6 for r in rows])
    w=np.array([float(r['Area']) for r in rows])
    require(0<len(rows)<=180000 and np.isfinite(angles).all() and np.isfinite(sizes).all() and np.isfinite(w).all() and min(sizes)>0 and min(w)>0,'Invalid input grain table.')
    w/=w.sum()
    policy=e['initial_grains']
    if policy!='measured_cohorts':
        sizes[:]=np.average(sizes,weights=w if policy=='equal_spheres_area_mean' else None)
        w[:]=1/len(w)
    return path,angles,sizes,w

def lattice_shift(c):
    return 30.0 if c['crystal']['system'] in ('hexagonal','trigonal') else 0.0

def crystal_text(c):
    x=c['crystal']
    a,b,cc=x['lattice_lengths_angstrom']
    mu=shear_modulus(c)/1e6
    nu=x['poisson_ratio']
    lam=2*mu*nu/(1-2*nu)
    C=np.eye(6)*mu
    C[:3,:3]=lam+np.eye(3)*2*mu
    modes=c['slip_modes']
    lines=[x['name'],x['system'].upper(),' '.join(f'{v:.15g}' for v in [1,b/a,cc/a,*x['lattice_angles_deg']]),
           'Auxiliary isotropic stiffness (MPa)',*[' '.join(f'{v:.12g}' for v in row) for row in C],
           'Isothermal expansion','0 0 0 0 0 0','Bidirectional slip modes',str(len(modes)),str(len(modes)),
           ' '.join(str(i+1) for i in range(len(modes)))]
    for i,m in enumerate(modes,1):
        lines.extend([m['name'],f'{i} {len(m["systems"])} 1 0',*[' '.join(map(str,row)) for row in m['systems']]])
    return '\n'.join(lines+['KM-DDRX','41','1','drx.nml'])+'\n'

def write_case(c,rate,case,stock_path=None):
    x,d,e=c['crystal'],c['drx'],c['experiment']
    path,angles,sizes,weights=grains(e)
    stock_material=None
    from .stock import uses_original
    if uses_original(c):
        from .stock import material_path, inspect_material
        require(len(weights)<=16000,'Original VPSC has a limit of 16,000 grains.')
        source=stock_path or material_path(c)
        shutil.copyfile(source,case/'magnesium.sx')
        stock_material=inspect_material(case/'magnesium.sx')
        shift=30.0 if stock_material['system'] in ('hexagonal','trigonal') else 0.0
    else:
        shift=lattice_shift(c)
    if e['orientation_convention']=='aztec':
        angles[:,2]=(angles[:,2]-shift)%360
    lines=['Bunge degrees, one phase, compression Z','','',f'B {len(weights)}']
    lines+=[' '.join(f'{v:.12g}' for v in [*a,w]) for a,w in zip(angles,weights)]
    (case/'initial.tex').write_text('\n'.join(lines)+'\n',encoding='ascii')
    de=e['max_strain_increment']
    if not stock_material:
        (case/'initial_sizes.dat').write_text(str(len(sizes))+'\n'+'\n'.join(f'{s:.15g}' for s in sizes)+'\n')
        if d['enabled']:
            surface_rate=d['nucleation_coefficient']*rate**(d['rate_exponent_q']-1)*math.pi*max(sizes)**2
            if surface_rate>0: de=min(de,.65/surface_rate)
    steps=max(1,math.ceil(e['end_strain']/de-1e-10))
    de=e['end_strain']/steps
    T=e['temperature_C']+273.15
    for name in ('VPSC8.in','compression.dat'):
        (case/name).write_text((TEMPLATES/name).read_text().format(steps=steps,increment=f'{de:.15g}',T=f'{T:.15g}',rate=f'{rate:.15g}',half=f'{rate/2:.15g}'))
    if stock_material:
        targets=c['plots']['strains']
        (case/'stock_snapshots.in').write_text(f'{rate:.17g} {len(targets)}\n'+' '.join(f'{t:.17g}' for t in targets)+'\n')
        save_json(case/'effective_inputs.json',dict(engine='original_vpsc',steps=steps,increment=de,rate_per_s=rate,
                  grain_sha256=sha(path),material_source=str(source),material_sha256=sha(case/'magnesium.sx'),material=stock_material))
        return steps
    modes=c['slip_modes']
    def values(key):
        return [m[key] for m in modes]
    config=dict(temperature=T,tau0=0,k1=modes[0]['storage_k1_per_m'],k2=modes[0]['recovery_k2'],
        rhoc=d['threshold_density_per_m2'],nuc_c=d['nucleation_coefficient'],mobility=d['mobility_m4_per_J_s'],
        burgers=x['burgers_reference_m'],chi=x['chi'],line_factor=d['line_energy_factor'],
        rho_initial=sum(values('initial_density_per_m2')),rho_nucleus=sum(values('nucleus_density_per_m2')),
        diameter_initial=float(np.mean(sizes**3)**(1/3)),diameter_nucleus=d['nucleus_diameter_m'],q=d['rate_exponent_q'],
        dt_limit=de/rate*1.00001,strain_rate=rate,seed_value=e['seed'],enabled=int(d['enabled']),
        snapshot_every=steps,slip_exponent=x['slip_exponent'],mu_reference=x['mu_reference_Pa'],mu_slope=x['mu_slope_per_K'],
        reference_temperature=x['mu_reference_temperature_K'],reference_macro_rate=x['recovery_reference_rate_per_s'],
        mode_tau=values('friction_MPa'),mode_b_ratio=[math.sqrt(1+(x['lattice_lengths_angstrom'][2]/x['lattice_lengths_angstrom'][0])**2)
            if m['burgers_ratio']=='hcp_c_plus_a' else m['burgers_ratio'] for m in modes],
        initial_size_file='initial_sizes.dat',mode_chi_ratio=values('chi_ratio'),mode_k1=values('storage_k1_per_m'),
        mode_k2=values('recovery_k2'),mode_recovery_rate_m=values('recovery_rate_exponent'),mode_rho_initial=values('initial_density_per_m2'),
        mode_rho_nucleus=values('nucleus_density_per_m2'),orientation_rule=int(d['orientation_rule']=='gaussian'),
        orientation_spread=d['orientation_spread_deg'],snapshot_strains=c['plots']['strains'])
    def fmt(v):
        return ', '.join(f'{i:.15g}' for i in v) if isinstance(v,list) else repr(v)
    (case/'drx.nml').write_text('&drx\n'+''.join(f' {k} = {fmt(v)}\n' for k,v in config.items())+'/\n')
    (case/'magnesium.sx').write_text(crystal_text(c),encoding='ascii')
    save_json(case/'effective_inputs.json',dict(steps=steps,increment=de,rate_per_s=rate,grain_sha256=sha(path),native_namelist=config))
    return steps


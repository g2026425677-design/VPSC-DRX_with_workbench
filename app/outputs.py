"""Output choices and lossless exports from completed per-rate solver histories."""
import csv
import json
from pathlib import Path
import shutil

from .model import require, save_json


FIGURES = {
    'stress': 'Axial stress–strain',
    'microstructure': 'DRX fraction, mean size and total density',
    'density': 'Dislocation density by family',
    'activity': 'Slip and twin activity',
    'pf': 'Pole figures (PF)',
    'ipf': 'Inverse pole figures (IPF)',
}
DATA = {
    **FIGURES,
    'grains': 'Grain states and orientations at saved strains',
    'cycles': 'DRX cycles and generation kinetics',
    'native': 'All native .OUT and solver/reporting CSV files',
}


def output_choices(settings):
    # A settings file from an earlier release produced all figures.
    options = settings.get('outputs', {'figures': list(FIGURES), 'data': []})
    require(isinstance(options, dict) and set(options) == {'figures', 'data'},
            'Outputs must contain figures and data lists.')
    for key, available in [('figures', FIGURES), ('data', DATA)]:
        values = options[key]
        require(isinstance(values, list) and all(isinstance(v, str) and v in available for v in values)
                and len(values) == len(set(values)), f'Invalid or duplicate {key} output selection.')
    return options


def subset_csv(source, target, fields):
    # Keep the original decimal strings, including all digits the solver wrote.
    with source.open(newline='') as src, target.open('w', newline='') as dst:
        reader = csv.DictReader(src, skipinitialspace=True)
        writer = csv.writer(dst)
        writer.writerow(fields)
        writer.writerows([row[field].strip() for field in fields] for row in reader)


def write_microstructure(case, dest, stock):
    with (case / 'MG_HISTORY.csv').open(newline='') as src, dest.open('w', newline='') as out:
        history = csv.DictReader(src, skipinitialspace=True)
        writer = csv.writer(out)
        if stock is None:
            with (case / 'DRX_HISTORY.csv').open(newline='') as src_drx:
                drx = csv.DictReader(src_drx, skipinitialspace=True)
                # The native strain and stress in this table are equivalent quantities.
                columns = [k for k in drx.fieldnames if k not in ('strain', 'stress_MPa')]
                writer.writerow(['axial_strain', *columns])
                writer.writerows([h['axial_strain'], *[d[k] for k in columns]]
                                 for h, d in zip(history, drx, strict=True))
        else:
            size = stock['hardening_grain_size_um']
            fields = ['axial_strain', 'drx_fraction']
            if size is not None:
                fields.append('hardening_grain_size_um')
            density = stock['hardening_law'] in (20, 23)
            if density:
                fields.append('rho_mean_m-2')
            writer.writerow(fields)
            def row(h):
                return [h['axial_strain'], 0] + ([size] if size is not None else [])
            if density:
                with (case / 'STOCK_DENSITY.csv').open(newline='') as src_rho:
                    for h, d in zip(history, csv.DictReader(src_rho, skipinitialspace=True), strict=True):
                        writer.writerow([*row(h), d['rho_total']])
            else:
                writer.writerows(row(h) for h in history)


def write_system_key(path, modes):
    with path.open('w', encoding='utf-8-sig', newline='') as dst:
        writer = csv.writer(dst)
        writer.writerow(['system', 'family', 'family_name', 'type', 'plane_indices', 'direction_indices'])
        index = 0
        for family, mode in enumerate(modes, 1):
            for system in mode['systems']:
                index += 1
                n = len(system) // 2
                writer.writerow([index, family, mode['name'], 'twin' if mode.get('twin') else 'slip',
                                 ' '.join(map(str, system[:n])), ' '.join(map(str, system[n:]))])


def export_rate(case, dest, selected, modes, stock):
    dest.mkdir(parents=True)
    notes = []
    if 'stress' in selected:
        subset_csv(case / 'MG_HISTORY.csv', dest / 'stress_strain.csv', ['axial_strain', 'axial_stress_MPa'])
    if 'microstructure' in selected:
        write_microstructure(case, dest / 'microstructure.csv', stock)
        if stock:
            notes.append('Original VPSC: DRX is off. Any hardening_grain_size_um is the material-file constant.')
    if 'density' in selected:
        if not stock or stock['hardening_law'] in (20, 23):
            source = 'STOCK_DENSITY.csv' if stock else 'MG_MODE_HISTORY.csv'
            shutil.copyfile(case / source, dest / 'density_by_family.csv')
    if 'activity' in selected:
        subset_csv(case / 'MG_HISTORY.csv', dest / 'family_activity.csv',
                   ['axial_strain', *[f'family_{i}' for i in range(1, len(modes)+1)]])
        shutil.copyfile(case / 'SYSTEM_ACTIVITY.csv', dest / 'system_activity.csv')
    if {'activity', 'density'} & selected:
        write_system_key(dest / 'system_key.csv', modes)
    if 'grains' in selected:
        grains = dest / 'grains'
        grains.mkdir()
        snapshots = case / 'MG_SNAPSHOTS.csv'
        shutil.copyfile(snapshots, grains / snapshots.name)
        with snapshots.open(newline='') as src:
            for row in csv.DictReader(src):
                filename = row['filename'].strip()
                require(Path(filename).name == filename and filename.startswith('MG_STATE_'),
                        'Invalid grain snapshot filename.')
                shutil.copyfile(case / filename, grains / filename)
    if 'cycles' in selected:
        if stock:
            notes.append('Original VPSC has no DRX generation history.')
        else:
            for name in ('DRX_CYCLES.csv', 'MG_CYCLE_KINETICS.csv'):
                shutil.copyfile(case / name, dest / name)
    if 'native' in selected:
        native = dest / 'native'
        native.mkdir()
        for path in sorted([*case.glob('*.OUT'), *case.glob('*.csv')]):
            shutil.copyfile(path, native / path.name)
    return notes


def export_outputs(c, run, dest):
    options = output_choices(c)
    figures, data = set(options['figures']), set(options['data'])
    notes = {}
    for rate in c['experiment']['rates_per_s']:
        case = run / f'rate_{rate}'
        meta = json.loads((case / 'effective_inputs.json').read_text(encoding='utf8'))
        stock = meta.get('material') if meta.get('engine') == 'original_vpsc' else None
        availability = []
        if data:
            modes = stock['modes'] if stock else c['slip_modes']
            availability = export_rate(case, dest / 'source_data' / case.name, data, modes, stock)
        if 'density' in (figures | data) and stock and stock['hardening_law'] not in (20, 23):
            availability.append('The selected original hardening law has no dislocation-density state.')
        notes[case.name] = availability
    if figures or {'pf', 'ipf'} & data:
        from .plots import render
        render(c, run, dest, figures, data)
    save_json(dest / 'output_manifest.json', dict(source_run=str(run.resolve()), selected=options,
              notes=notes, files=sorted(p.relative_to(dest).as_posix() for folder in ('figures', 'source_data')
                                        for p in (dest / folder).rglob('*') if p.is_file())))

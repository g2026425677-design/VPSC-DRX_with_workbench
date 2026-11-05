"""Plot native histories and official MTEX PF/IPF densities using current grain volumes."""
import csv
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
import numpy as np
from .model import ROOT, SYMMETRY, lattice_shift, require, save_json

if sys.platform=='win32':
    sys.path.insert(0,str(ROOT/'dependencies'))
import mtex as mt

def read(path):
    return np.genfromtxt(path,delimiter=',',names=True,deletechars='',encoding='ascii',ndmin=1)

def chart(ax,ylabel):
    ax.set_xlabel('Axial true compressive strain')
    ax.set_ylabel(ylabel)
    ax.grid(alpha=.2)
    ax.spines[['top','right']].set_visible(False)

def histories(c,out,dest,stock,figures):
    rates=c['experiment']['rates_per_s']
    modes=c['slip_modes']
    if 'stress' in figures:
        stress,sa=plt.subplots(figsize=(7.3,4.6),layout='constrained')
    if 'microstructure' in figures:
        micro,axes=plt.subplots(2,2,figsize=(11.4,7.4),layout='constrained')
    for rate,color in zip(rates,plt.cm.viridis(np.linspace(.12,.83,len(rates)))):
        case=out/f'rate_{rate}'
        h=read(case/'MG_HISTORY.csv')
        rho_path=case/('STOCK_DENSITY.csv' if stock else 'MG_MODE_HISTORY.csv')
        rho=read(rho_path) if figures & {'microstructure','density'} and rho_path.exists() else None
        strain=h['axial_strain']
        label=f'{rate:g} s⁻¹'
        if 'stress' in figures:
            sa.plot(strain,h['axial_stress_MPa'],label=label,color=color)
        if 'microstructure' in figures:
            if stock:
                axes.flat[0].plot(strain,np.zeros_like(strain),label=label,color=color)
                if stock['hardening_grain_size_um'] is not None:
                    axes.flat[1].plot(strain,np.full_like(strain,stock['hardening_grain_size_um']),label=label,color=color)
                if rho is not None:
                    axes.flat[2].plot(strain,rho['rho_total'],label=label,color=color)
            else:
                d=read(case/'DRX_HISTORY.csv')
                for ax,key in zip(axes.flat,['drx_fraction','mean_diameter_volume_um','rho_mean_m-2','mean_diameter_number_um']):
                    ax.plot(strain,d[key],label=label,color=color)
        if 'density' in figures and rho is not None:
            fig,ax=plt.subplots(figsize=(7,4.4),layout='constrained')
            for i,m in enumerate(modes,1):
                if f'rho_{i}' in rho.dtype.names:
                    ax.plot(strain,rho[f'rho_{i}'],label=m['name'])
            if stock and stock['hardening_law']==20:
                ax.plot(strain,rho['rho_debris'],label='Debris')
            ax.plot(strain,rho['rho_total'],color='#262626',ls='--',label='Total')
            chart(ax,'Volume-mean dislocation density (m⁻²)')
            ax.legend(fontsize=8)
            ax.set_title(f'Mode-resolved density  |  {label}')
            fig.savefig(dest/f'dislocation_modes_rate_{rate}.png',dpi=160)
            plt.close(fig)
        if 'activity' not in figures:
            continue
        act=read(case/'SYSTEM_ACTIVITY.csv')
        fig,panels=plt.subplots(len(modes)+1,1,figsize=(9.4,2.45*(len(modes)+1)),layout='constrained',squeeze=False)
        fas=panels[0,0]
        first=1
        for i,(m,ax) in enumerate(zip(modes,panels[1:,0]),1):
            fas.plot(strain,100*h[f'family_{i}'],label=m['name'])
            count=len(m['systems'])
            colors=plt.cm.turbo(np.linspace(.04,.93,count))
            for j,(row,sc) in enumerate(zip(m['systems'],colors),first):
                half=len(row)//2
                fmt=lambda v:' '.join(map(str,v))
                ax.plot(strain,100*act[f'system_{j}'],color=sc,label=f'S{j}   ({fmt(row[:half])}) [{fmt(row[half:])}]')
            first+=count
            chart(ax,'Share of total shear (%)' if stock else 'Share of total slip (%)')
            ax.set_title(m['name'],loc='left')
            ax.legend(fontsize=7,loc='upper left',bbox_to_anchor=(1.01,1),ncol=2 if count>12 else 1)
        chart(fas,'Share of total shear (%)' if stock else 'Share of total slip (%)')
        fas.legend(fontsize=8)
        fas.set_title(('Original VPSC slip and twin activity' if stock else 'Slip-family and individual-system activity')+f'  |  {label}',loc='left')
        fig.savefig(dest/f'slip_activity_rate_{rate}.png',dpi=150)
        plt.close(fig)
    heading=f"{c['experiment']['temperature_C']:g} °C, uniaxial Z compression"
    if stock:
        heading=f'Original VPSC · law {stock["hardening_law"]}  |  '+heading
        if 'microstructure' in figures:
            axes.flat[0].set_ylim(0,1)
    if 'stress' in figures:
        chart(sa,'Axial true compressive stress (MPa)')
        sa.set_title(heading)
        sa.legend()
        stress.savefig(dest/'stress_strain.png',dpi=180)
        plt.close(stress)
    if 'microstructure' not in figures:
        return
    labels=['Historical DRX volume fraction','Hardening grain size from .sx (µm)' if stock else 'Volume-weighted diameter (µm)','Volume-mean total density (m⁻²)','Number-weighted diameter (µm)']
    for i,(ax,label) in enumerate(zip(axes.flat,labels)):
        if stock and not ax.lines:
            ax.axis('off')
            note={1:'No physical diameter reported for this law.',2:'This material law has no density state.',3:'Original VPSC\nNo evolving DRX grain population.'}[i]
            ax.text(.5,.5,note,transform=ax.transAxes,ha='center',va='center',color='#555')
        else:
            chart(ax,label)
            ax.legend(fontsize=8)
    micro.suptitle(heading)
    micro.savefig(dest/'microstructure.png',dpi=180)
    plt.close(micro)

def field(v,w,kernel):
    density=mt.calcDensity(v,weights=w,kernel=kernel,antipodal=True,bandwidth=kernel.bandwidth)
    require(abs(float(density.mean())-1)<1e-7,'MTEX density normalization failed.')
    return density

def projection_options(kind):
    return dict(upper=True, antipodal=True, projection='earea' if kind=='PF' else 'eangle',
                how2plot=mt.plottingConvention('y↑→x'), fundamentalRegion=kind=='IPF')


def export_texture_grid(density, kind, path):
    from mtex.plotting import makeSphericalProjection
    projection, = makeSphericalProjection(density, **projection_options(kind))
    vectors, x, y = projection.grid(mt.degree)
    directions = mt.vector3d(vectors.data, frame=density.frame, isNormalized=True)
    mrd = np.real(np.asarray(density(directions))).reshape(x.shape)
    valid = np.isfinite(x) & np.isfinite(y) & np.isfinite(mrd)
    require(valid.any(), 'Empty MTEX density grid.')
    grid = np.column_stack([x[valid], y[valid], vectors.data[valid], mrd[valid]])
    np.savetxt(path, grid, delimiter=',', fmt='%.17g', comments='',
               header='projection_x,projection_y,direction_x,direction_y,direction_z,MRD')


def texture(c,out,dest,figure_keys,data_keys):
    figure_dir = dest / 'figures'
    kinds = figure_keys | data_keys
    mt.setMTEXpref('threads',1)
    x=c['crystal']
    p=c['plots']
    cs=mt.crystalFrame(SYMMETRY[x['system']],x['lattice_lengths_angstrom'],np.asarray(x['lattice_angles_deg'])*mt.degree,
                       mineral=x['name'])
    kernel=mt.S2DeLaValleePoussinKernel(halfwidth=p['halfwidth_deg']*mt.degree)
    records=[]
    maps=[]
    pfamilies=[]
    planes = p['pole_planes'] if 'pf' in kinds else []
    directions = p['ipf_specimen_directions'] if 'ipf' in kinds else []
    for plane in planes:
        idx='hkil' if len(plane)==4 else 'hkl'
        pole=mt.Miller(*plane,cs,idx).normalize()
        pfamilies.append(pole.symmetrise(unique=True,antipodal=True))
    for rate in c['experiment']['rates_per_s']:
        case=out/f'rate_{rate}'
        with (case/'MG_SNAPSHOTS.csv').open() as f:
            snaps=list(csv.DictReader(f))
        grid_rows = []
        if data_keys:
            grid_dir = dest / 'source_data' / case.name / 'texture'
            grid_dir.mkdir(parents=True, exist_ok=True)
        for j,s in enumerate(snaps,1):
            state=read(case/s['filename'].strip())
            angles=np.column_stack([state[k] for k in ('phi1','Phi','phi2')])
            if c['experiment']['orientation_convention']!='archived':
                angles[:,2]=(angles[:,2]+lattice_shift(c))%360
            ori=mt.orientation.byEuler(*(angles*mt.degree).T,cs)
            w=state['volume_fraction']
            require(np.all(w>0) and abs(w.sum()-1)<1e-8,'Invalid native snapshot volume weights.')
            densities=[]
            for plane,family in zip(planes,pfamilies):
                vectors=np.concatenate([(ori*h).normalize().data for h in family])
                dens=field(mt.vector3d(vectors),np.tile(w,len(family)),kernel)
                densities.append(('PF',plane,dens))
            for direction in directions:
                v=mt.vector3d(np.asarray(direction,float)).normalize()
                dens=field(ori.inv()*v,w,kernel)
                densities.append(('IPF',direction,dens))
            visible = []
            for k,(kind,indices,dens) in enumerate(densities,1):
                peak,_=mt.max(dens,resolution=1*mt.degree,accuracy=.05*mt.degree)
                rec=dict(kind=kind,indices=indices,rate=rate,target=float(s['target_strain']),actual=float(s['actual_axial_strain']),
                         grains=len(w),peak_MRD=float(peak),snapshot=s['filename'].strip())
                if kind.lower() in data_keys:
                    grid_name = f'{kind}_state_{j:02d}_map_{k:02d}.csv'
                    export_texture_grid(dens, kind, grid_dir / grid_name)
                    grid_rows.append(dict(**rec,indices_text=' '.join(map(str,indices)),file=grid_name,
                                          halfwidth_deg=p['halfwidth_deg'],resolution_deg=1,MTEX_version=mt.__version__))
                if kind.lower() in figure_keys:
                    visible.append((rec, (kind, indices, dens)))
                    records.append(rec)
            if visible:
                recs, shown = zip(*visible)
                maps.append((rate,j,recs,shown))
        if grid_rows:
            columns = ['file','kind','indices_text','rate','target','actual','peak_MRD','grains','snapshot',
                       'halfwidth_deg','resolution_deg','MTEX_version']
            with (grid_dir/'texture_index.csv').open('w',newline='') as target:
                writer=csv.DictWriter(target,fieldnames=columns,extrasaction='ignore')
                writer.writeheader()
                writer.writerows(grid_rows)
    if not maps:
        return
    scale=max(1,np.ceil(max(r['peak_MRD'] for r in records)))
    for rate,j,recs,densities in maps:
        fig,axes=plt.subplots(1,len(densities),figsize=(4.1*len(densities),4.75),squeeze=False)
        for ax,rec,(kind,indices,density) in zip(axes[0],recs,densities):
            mt.plot(density, parent=ax, **projection_options(kind), pcolor=True, colorRange=[0,scale],
                    colormap='viridis', resolution=mt.degree)
            text=' '.join(f'{v:g}' for v in indices)
            ax.set_title(f'{kind}  {text}\nMaximum {rec["peak_MRD"]:.2f} MRD',fontsize=11)
        fig.subplots_adjust(top=.78,bottom=.22,left=.055,right=.95,wspace=.24)
        fig.suptitle(f"{c['experiment']['temperature_C']:g} °C  |  {rate:g} s⁻¹  |  axial strain {recs[0]['actual']:.5g}",fontsize=14)
        cb=fig.colorbar(ScalarMappable(norm=Normalize(0,scale),cmap='viridis'),cax=fig.add_axes([.23,.12,.54,.022]),orientation='horizontal')
        cb.set_label(f"Pole density (MRD)   ·   MTEX halfwidth {p['halfwidth_deg']:g}°")
        fig.savefig(figure_dir/f'texture_rate_{rate}_state_{j:02d}.png',dpi=175,bbox_inches='tight',pad_inches=.16)
        plt.close(fig)
    save_json(figure_dir/'texture_metrics.json',dict(MTEX_version=mt.__version__,shared_color_max_MRD=float(scale),records=records))

def render(c,out,dest,figures,data):
    figure_dir=dest/'figures'
    if figures:
        figure_dir.mkdir(exist_ok=True)
    first=c['experiment']['rates_per_s'][0]
    meta=json.loads((out/f'rate_{first}'/'effective_inputs.json').read_text(encoding='utf8'))
    stock=meta.get('material') if meta.get('engine')=='original_vpsc' else None
    if stock:
        # Use the geometry, mode order and twin IDs from the .sx used by VPSC.
        c=dict(c,crystal=stock,slip_modes=stock['modes'])
    texture_keys = {'pf', 'ipf'}
    if figures - texture_keys:
        histories(c,out,figure_dir,stock,figures)
    if (figures | data) & texture_keys:
        texture(c,out,dest,figures & texture_keys,data & texture_keys)


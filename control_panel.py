"""Desktop parameter editor and launcher for the calibrated VPSC-DDRX solver."""
import copy
import json
import os
import subprocess
import sys
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from datetime import datetime
from pathlib import Path
from app.model import ROOT, validate, grains
from app.outputs import FIGURES, DATA, output_choices
from app.stock import uses_original, material_path, inspect_material, DEFAULT_SX

def as_text(v):
    return json.dumps(v,ensure_ascii=False) if isinstance(v,(list,bool)) else str(v)

def parse(text, sample):
    text=text.strip()
    if isinstance(sample,str):
        return text
    return json.loads(text)

class ModeEditor:
    def __init__(self,frame,modes):
        self.modes=copy.deepcopy(modes)
        self.index=0
        self.fields={}
        bar=ttk.Frame(frame)
        bar.grid(row=1,column=0,columnspan=2,sticky='w',pady=5)
        self.menu=ttk.Combobox(bar,state='readonly',width=34)
        self.menu.pack(side='left')
        self.menu.bind('<<ComboboxSelected>>',self.switch)
        ttk.Button(bar,text='Add family',command=self.add).pack(side='left',padx=5)
        ttk.Button(bar,text='Remove family',command=self.remove).pack(side='left')
        row=2
        for key in self.modes[0]:
            if key=='systems':
                continue
            ttk.Label(frame,text=key.replace('_',' ')).grid(row=row,column=0,sticky='w',pady=2,padx=(0,10))
            self.fields[key]=tk.StringVar()
            ttk.Entry(frame,textvariable=self.fields[key]).grid(row=row,column=1,sticky='ew',pady=2)
            row+=1
        ttk.Label(frame,text='Slip systems · one plane and direction per line').grid(row=row,column=0,columnspan=2,sticky='w',pady=5)
        self.rows=tk.Text(frame,height=6,wrap='none',font=('Consolas',10))
        self.rows.grid(row=row+1,column=0,columnspan=2,sticky='nsew')
        scrollbar=ttk.Scrollbar(frame,command=self.rows.yview)
        scrollbar.grid(row=row+1,column=2,sticky='ns')
        self.rows.configure(yscrollcommand=scrollbar.set)
        frame.rowconfigure(row+1,weight=1)
        self.display(0)

    def display(self,index):
        self.index=index
        m=self.modes[index]
        self.menu.configure(values=[f'{i+1}. {v["name"]}' for i,v in enumerate(self.modes)])
        self.menu.current(index)
        for k,var in self.fields.items():
            var.set(as_text(m[k]))
        self.rows.delete('1.0','end')
        self.rows.insert('1.0','\n'.join('   '.join(map(str,v)) for v in m['systems']))

    def store(self):
        m={}
        for k,var in self.fields.items():
            value=var.get().strip()
            m[k]=value if k=='name' or (k=='burgers_ratio' and value=='hcp_c_plus_a') else float(value)
        m['systems']=[list(map(int,line.split())) for line in self.rows.get('1.0','end').splitlines() if line.strip()]
        self.modes[self.index]=m
        return self.modes

    def switch(self,event):
        selected=self.menu.current()
        try:
            self.store()
            self.display(selected)
        except ValueError as exc:
            self.menu.current(self.index)
            messagebox.showerror('Invalid slip family',str(exc))

    def add(self):
        try:
            self.store()
            if len(self.modes)>=8:
                raise ValueError('The native limit is eight slip families.')
            fresh=copy.deepcopy(self.modes[self.index])
            fresh.update(name='New family',systems=[])
            self.modes.append(fresh)
            self.display(len(self.modes)-1)
        except ValueError as exc:
            messagebox.showerror('Cannot add family',str(exc))

    def remove(self):
        if len(self.modes)>1:
            self.modes.pop(self.index)
            self.display(min(self.index,len(self.modes)-1))

class Panel:
    def __init__(self, window):
        self.window=window
        self.path=ROOT/'settings.json'
        self.proc=None
        self.output=None
        self.log=None
        self.fields={}
        self.window.title('VPSC–DRX  |  Compression workbench')
        self.window.geometry('1010x785')
        self.window.minsize(770,640)
        head=ttk.Frame(window,padding=(18,14,18,6))
        head.pack(fill='x')
        ttk.Label(head,text='VPSC–DRX compression',font=('Segoe UI',18,'bold')).pack(anchor='w')
        ttk.Label(head,text='Edit a material, select the strain states, then run the native solver and MTEX.').pack(anchor='w')
        top=ttk.Frame(window,padding=(18,2))
        top.pack(fill='x')
        for label,cmd in [('Load settings…',self.load),('Save as…',self.save_as),('Restore calibrated values',self.reset)]:
            ttk.Button(top,text=label,command=cmd).pack(side='left',padx=(0,7))
        self.notebook=ttk.Notebook(window)
        self.notebook.pack(fill='both',expand=True,padx=18,pady=8)
        bottom=ttk.Frame(window,padding=(18,6,18,12))
        bottom.pack(fill='x')
        self.run_button=ttk.Button(bottom,text='Save & run simulation',command=self.run)
        self.run_button.pack(side='left')
        ttk.Button(bottom,text='Open results',command=self.open_results).pack(side='left',padx=8)
        ttk.Button(bottom,text='Open guide',command=lambda: os.startfile(ROOT/'README.md')).pack(side='left')
        self.status=tk.StringVar()
        ttk.Label(bottom,textvariable=self.status).pack(side='left',padx=10)
        self.show(json.loads(self.path.read_text(encoding='utf-8')))
        self.window.protocol('WM_DELETE_WINDOW',self.close)

    def show(self,c):
        self.config=copy.deepcopy(c)
        self.fields={}
        for tab in self.notebook.tabs():
            self.notebook.nametowidget(tab).destroy()
        tabs=[('experiment','Test & grains','Positive true compression along Z. Lists use [comma-separated values].'),
              ('crystal','Crystal','Lattice lengths in Å, angles in degrees. The hot shear modulus is in Pa.'),
              ('slip_modes','Slip & hardening','Each row is plane indices followed by direction indices. Change a family’s parameters below.'),
              ('drx','Recrystallization','DRX on uses the calibrated paper K–M law. DRX off uses the hardening law selected below. Native-law coefficients are edited in Original material.'),
              ('plots','Texture settings','Enter up to 16 ordered snapshot strains. Directions for inverse pole figures use specimen X, Y, Z.')]
        for section,label,note in tabs:
            frame=ttk.Frame(self.notebook,padding=14)
            self.notebook.add(frame,text=label)
            frame.columnconfigure(1,weight=1)
            ttk.Label(frame,text=note,wraplength=840).grid(row=0,column=0,columnspan=2,sticky='w',pady=(0,8))
            if section=='slip_modes':
                self.mode_editor=ModeEditor(frame,c[section])
                continue
            choices={'system':['hexagonal','cubic','tetragonal','orthorhombic','trigonal'],
                     'orientation_convention':['archived','aztec','vpsc'],
                     'initial_grains':['measured_cohorts','equal_spheres_number_mean','equal_spheres_area_mean'],
                     'orientation_rule':['gaussian','random'],'enabled':['true','false'],
                     'hardening_when_off':['original_dd_rev','original_voce','original_material','paper_km']}
            values=dict(c[section])
            if section=='drx':
                values.setdefault('original_material_file',DEFAULT_SX)
                values.setdefault('hardening_when_off','original_dd_rev')
            for row,(key,v) in enumerate(values.items(),1):
                ttk.Label(frame,text=key.replace('_',' ')).grid(row=row,column=0,sticky='w',pady=4,padx=(0,14))
                var=tk.StringVar(value=as_text(v))
                self.fields[section,key]=var
                if key in choices:
                    entry=ttk.Combobox(frame,textvariable=var,values=choices[key],state='readonly')
                    if key=='hardening_when_off':
                        entry.bind('<<ComboboxSelected>>',self.select_native_law)
                else:
                    entry=ttk.Entry(frame,textvariable=var)
                entry.grid(row=row,column=1,sticky='ew',pady=3)
        native=ttk.Frame(self.notebook,padding=14)
        self.notebook.add(native,text='Original material')
        ttk.Label(native,text='Edit the original .sx coefficients and slip/twin definitions here. Save & run saves edits as a separate material file. The selected law number and its required parameter block must agree.',wraplength=860).pack(anchor='w')
        bar=ttk.Frame(native)
        bar.pack(fill='x',pady=8)
        for label,command in [('Browse .sx…',self.browse_native),('Reload selected file',self.load_native)]:
            ttk.Button(bar,text=label,command=command).pack(side='left',padx=(0,7))
        self.native_info=tk.StringVar()
        ttk.Label(bar,textvariable=self.native_info).pack(side='left')
        body=ttk.Frame(native)
        body.pack(fill='both',expand=True)
        self.native_text=tk.Text(body,wrap='none',font=('Consolas',10),undo=True)
        scroll=ttk.Scrollbar(body,command=self.native_text.yview)
        scroll.pack(side='right',fill='y')
        self.native_text.pack(side='left',fill='both',expand=True)
        self.native_text.configure(yscrollcommand=scroll.set)
        self.build_outputs(c)
        self.load_native()
        for key in ('enabled','hardening_when_off'):
            self.fields['drx',key].trace_add('write',lambda *_: self.update_law_tabs())
        self.update_law_tabs()
        self.status.set(self.path.name)

    def build_outputs(self, c):
        frame = ttk.Frame(self.notebook, padding=14)
        self.notebook.add(frame, text='Outputs')
        options = output_choices(c)
        self.output_vars = {'figures': {}, 'data': {}}
        for col, title in enumerate(('Quantity', 'Figure (PNG)', 'Source data')):
            ttk.Label(frame, text=title).grid(row=0, column=col, sticky='w', padx=8, pady=(0,8))
        for row, (key, label) in enumerate(DATA.items(), 1):
            ttk.Label(frame, text=label).grid(row=row, column=0, sticky='w', padx=8, pady=3)
            for col, group in enumerate(('figures', 'data'), 1):
                if group == 'figures' and key not in FIGURES:
                    continue
                var = tk.BooleanVar(value=key in options[group])
                self.output_vars[group][key] = var
                ttk.Checkbutton(frame, variable=var).grid(row=row, column=col)
        row = len(DATA) + 1
        ttk.Label(frame, text='Source data are grouped by strain rate. PF/IPF data are MTEX density grids in MRD.',
                  wraplength=680).grid(row=row,column=0,columnspan=3,sticky='w',padx=8,pady=(12,4))
        self.export_button = ttk.Button(frame,text='Export from a completed run…',command=self.export_saved)
        self.export_button.grid(row=row+1,column=0,sticky='w',padx=8,pady=4)
        if self.proc and self.proc.poll() is None:
            self.export_button.state(['disabled'])
        ttk.Label(frame,text='Uses the selected run’s saved strain states and texture settings.',
                  wraplength=680).grid(row=row+2,column=0,columnspan=3,sticky='w',padx=8)

    def selected_outputs(self):
        return {group: [key for key,var in fields.items() if var.get()]
                for group, fields in self.output_vars.items()}

    def export_saved(self):
        name = filedialog.askdirectory(title='Select a run containing settings_used.json and rate folders',
                                       initialdir=ROOT/'results')
        if not name:
            return
        try:
            run = Path(name)
            if not (run/'settings_used.json').is_file():
                raise ValueError('Select the run folder containing settings_used.json.')
            stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            base = run/'exports'
            base.mkdir(exist_ok=True)
            options = self.selected_outputs()
            self.launch(['--export-from',str(run),'--figures',*options['figures'],'--data',*options['data']],
                        base/stamp, 'Exporting saved simulation…')
        except Exception as exc:
            messagebox.showerror('Cannot export',str(exc))

    def update_law_tabs(self):
        stock=self.fields['drx','enabled'].get()=='false' and self.fields['drx','hardening_when_off'].get()!='paper_km'
        for index in (1,2):
            self.notebook.tab(index,state='disabled' if stock else 'normal')

    def load_native(self):
        path=(ROOT/Path(self.fields['drx','original_material_file'].get())).resolve()
        try:
            info=inspect_material(path)
            content=path.read_bytes().decode('latin1')
            self.native_text.delete('1.0','end')
            self.native_text.insert('1.0',content.replace('\r\n','\n').rstrip('\x1a'))
            self.native_loaded_text=self.native_text.get('1.0','end-1c')
            self.native_loaded_path=path
            self.native_info.set(f'IHARDLAW {info["hardening_law"]} · {path.name}')
        except Exception as exc:
            self.native_info.set(f'Cannot load material · {exc}')
            self.native_loaded_path=None

    def browse_native(self):
        value=filedialog.askopenfilename(initialdir=ROOT/'data/original_vpsc',filetypes=[('VPSC material','*.sx'),('All files','*')])
        if value:
            self.fields['drx','original_material_file'].set(value)
            self.fields['drx','hardening_when_off'].set('original_material')
            self.load_native()

    def select_native_law(self,event=None):
        law=self.fields['drx','hardening_when_off'].get()
        presets={'original_dd_rev':'data/original_vpsc/Mg.sx','original_voce':'data/original_vpsc/Mg_voce.sx'}
        if law in presets:
            self.fields['drx','original_material_file'].set(presets[law])
            self.load_native()

    def store_native_edits(self):
        c=dict(self.config,drx=dict(self.config['drx'],enabled=self.fields['drx','enabled'].get()=='true',
                                   hardening_when_off=self.fields['drx','hardening_when_off'].get()))
        if not uses_original(c):
            return
        selected=(ROOT/Path(self.fields['drx','original_material_file'].get())).resolve()
        content=self.native_text.get('1.0','end-1c')
        if selected!=self.native_loaded_path:
            if content!=getattr(self,'native_loaded_text',''):
                raise ValueError('Reload the selected material path before editing its contents.')
            self.load_native()
            return
        if content==self.native_loaded_text:
            return
        dest=ROOT/'data/user_materials'
        dest.mkdir(exist_ok=True)
        path=dest/f'material_{datetime.now().strftime("%Y%m%d_%H%M%S_%f")}.sx'
        path.write_bytes((content.rstrip('\n')+'\n').encode('latin1'))
        inspect_material(path)
        self.fields['drx','original_material_file'].set(str(path.relative_to(ROOT)))
        self.load_native()

    def collect(self):
        c=copy.deepcopy(self.config)
        for (s,k),v in self.fields.items():
            try:
                c[s][k]=parse(v.get(),c[s].get(k,''))
            except ValueError as exc:
                raise ValueError(f'{k}: enter a number or JSON list.') from exc
        if not uses_original(c):
            c['slip_modes']=copy.deepcopy(self.mode_editor.store())
        c['outputs']=self.selected_outputs()
        validate(c)
        grains(c['experiment'])
        return c

    def save(self,path):
        self.store_native_edits()
        c=self.collect()
        tmp=path.with_suffix(path.suffix+'.tmp')
        tmp.write_text(json.dumps(c,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
        tmp.replace(path)
        self.config=c
        self.path=path
        self.status.set(f'Saved {path.name}')

    def load(self):
        p=filedialog.askopenfilename(initialdir=ROOT,filetypes=[('Settings','*.json')])
        if p:
            try:
                c=json.loads(Path(p).read_text(encoding='utf-8-sig'))
                validate(c)
                self.path=Path(p)
                self.show(c)
            except Exception as exc:
                messagebox.showerror('Invalid settings',str(exc))

    def save_as(self):
        p=filedialog.asksaveasfilename(initialdir=ROOT,defaultextension='.json')
        if p:
            try:
                self.save(Path(p))
            except Exception as exc:
                messagebox.showerror('Cannot save',str(exc))

    def reset(self):
        self.show(json.loads((ROOT/'app/default_settings.json').read_text(encoding='utf8')))
        self.status.set('Calibrated defaults loaded into editor')

    def run(self):
        try:
            self.save(self.path)
            stamp=datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            base=ROOT/'results'
            base.mkdir(exist_ok=True)
            self.launch(['--settings',str(self.path)],base/stamp,'Native calculation running…')
        except Exception as exc:
            messagebox.showerror('Cannot start simulation',str(exc))

    def launch(self, args, output, status):
        self.output=output
        self.log=output.parent/(output.name+'_launcher.log')
        with self.log.open('w',encoding='utf8') as f:
            self.proc=subprocess.Popen([sys.executable,'-X','utf8',str(ROOT/'simulate.py'),*args,'--output',str(output)],
                stdout=f,stderr=subprocess.STDOUT,cwd=ROOT,creationflags=getattr(subprocess,'CREATE_NO_WINDOW',0))
        self.run_button.state(['disabled'])
        self.export_button.state(['disabled'])
        self.status.set(status)
        self.window.after(1000,self.poll)

    def poll(self):
        code=self.proc.poll()
        if code is None:
            lines=self.log.read_text(encoding='utf8',errors='replace').splitlines()
            if lines and len(lines[-1])<95:
                self.status.set(lines[-1])
            self.window.after(2000,self.poll)
            return
        self.run_button.state(['!disabled'])
        self.export_button.state(['!disabled'])
        self.status.set('Completed. Use Open results to view the outputs.' if code==0 else 'Run failed. Open results or the launcher log.')
        if code:
            messagebox.showerror('Run failed',self.log.read_text(encoding='utf8',errors='replace')[-1400:])

    def open_results(self):
        folder=self.output or ROOT/'results/calibrated_example'
        if folder.is_dir():
            os.startfile(folder)

    def close(self):
        if self.proc and self.proc.poll() is None:
            if not messagebox.askyesno('Close editor?','The simulation will continue. Its results will be written to the selected folder.'):
                return
        self.window.destroy()

if __name__=='__main__':
    win=tk.Tk()
    Panel(win)
    win.mainloop()


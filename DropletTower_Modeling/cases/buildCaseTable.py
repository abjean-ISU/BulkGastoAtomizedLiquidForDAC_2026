#builds BGAL_cases.xlsx from BGAL Settings.xlsx and the ILCompare screening inputs
#inputs are typed values; derived columns are Excel formulas (same calculations as the BGAL Settings sheets);
#the Aspen parameter input columns are filled by runAll.m
#usage: python buildCaseTable.py "BGAL Settings.xlsx" BGAL_cases.xlsx [--overwrite]
import os, sys, openpyxl
from openpyxl.utils import get_column_letter
from openpyxl.styles import Font, PatternFill

src, out = sys.argv[1], sys.argv[2]
if os.path.exists(out) and '--overwrite' not in sys.argv:
    sys.exit(f'{out} exists and may hold MATLAB results; pass --overwrite to rebuild it')
wb_in = openpyxl.load_workbook(src, data_only=True, read_only=True)

ABS = {'Height':'towerHeight','Diameter':'towerDiameter','Orifice Density':'orificeDensity','Orifice Diameter':'d_orifice_mm',
       'Flowrate per Orifice':'Q_orifice','Temp':'T_C','Pressure':'P','ppmCO2':'ppm0_CO2','ppmO2':'ppm0_O2','ppmN2':'ppm0_N2',
       'gas velocity - SET':'Ug','Vapor Pressure':'P_sat_sorbent',
       'CO2 Capture eff':'eff_CO2_sheet','O2 Capture eff':'eff_O2_sheet','N2 Capture eff':'eff_N2_sheet'}
DES = dict(ABS, **{'FCO2 - Aspen Results':'F_CO2','FO2 - Aspen Results':'F_O2','FN2 - Aspen Results':'F_N2',
                   'CO2 desorption eff':'eff_CO2_sheet','O2 desorption eff':'eff_O2_sheet','N2 desorption eff':'eff_N2_sheet'})

def grab(ws, col, labels):
    vals = {}
    for r in ws.iter_rows(values_only=True):
        if len(r) > col+1 and isinstance(r[col], str):
            key = r[col].strip()
            if key in labels and labels[key] not in vals:
                vals[labels[key]] = r[col+1]
    return vals

#column name, units, type: 1 = typed input, 0 = formula, 2 = MATLAB output
cols = [('caseName','',1),('sorbent','',1),('stage','abs/des',1),
        ('towerHeight','m',1),('towerDiameter','m',1),('A_c','m2',0),('manifoldDiameter','m',0),
        ('orificeDensity','orifice/m2',1),('N_o','orifice',0),('d_orifice_mm','mm',1),('Q_orifice','m3/s',1),
        ('Q','m3/s',0),('v0','m/s',0),('T_C','C',1),('P','atm',1),('P_baseline','atm',1),('Ug_baseline','m/s',1),
        ('ppm0_CO2','ppm',1),('ppm0_O2','ppm',1),('ppm0_N2','ppm',1),
        ('F_CO2','kg/s',1),('F_O2','kg/s',1),('F_N2','kg/s',1),('c0_CO2','g/L',0),('c0_O2','g/L',0),('c0_N2','g/L',0),
        ('P_sat_sorbent','atm',1),('H_scale','-',1),('Ncol','-',1),
        ('eff_CO2_sheet','%',1),('eff_O2_sheet','%',1),('eff_N2_sheet','%',1),
        #Aspen parameter inputs
        ('Ug','m/s',2),('mdot_gas','kg/s',2),('mdot_liquid','kg/s',2),('P_flash','atm',0),('P_pump','atm',0),
        ('T_flash','C',0),('eff_CO2','%',2),('eff_O2','%',2),('eff_N2','%',2),('P_vapor','atm',0),('approachSat','%',0)]
L = {name: get_column_letter(i+1) for i, (name, _, _) in enumerate(cols)}

def formulas(r, stage, sorbent):
    c = lambda name: f'{L[name]}{r}'
    dP_pump = 1.55 if sorbent == 'water' else 1.65 #pump head above column P (atm), as in the BGAL Settings sheets
    f = {'A_c': f'=PI()*({c("towerDiameter")}/2)^2',
         'manifoldDiameter': f'={c("towerDiameter")}-0.1',
         'N_o': f'=CEILING({c("orificeDensity")}*PI()*({c("manifoldDiameter")}/2)^2,1)',
         'Q': f'={c("Q_orifice")}*{c("N_o")}',
         'v0': f'={c("Q")}/(PI()*({c("d_orifice_mm")}/1000/2)^2*{c("N_o")})',
         'P_flash': f'={c("P")}', 'T_flash': f'={c("T_C")}', 'P_vapor': f'={c("P_sat_sorbent")}'}
    if stage == 'abs':
        f.update(c0_CO2=0, c0_O2=0, c0_N2=0, P_pump=f'={c("P")}+{dP_pump}')
    else:
        f.update(c0_CO2=f'={c("F_CO2")}/{c("Q")}', c0_O2=f'={c("F_O2")}/{c("Q")}', c0_N2=f'={c("F_N2")}/{c("Q")}',
                 approachSat=f'=IF(COUNT({c("eff_CO2")}:{c("eff_N2")})=3,AVERAGE({c("eff_CO2")}:{c("eff_N2")}),"")')
    return f

pressureCases = {'W3atm','W6atm','W15atm','W20atm','IL3atm','IL6atm','IL15atm','IL20atm'}
H_scale = {'IL High Abs':0.5, 'IL Low Abs':1.5}
rows = []

#IL screening (inputs from the *_ILCompare.m scripts; 10,000 orifices/m2 on the 11.9 m manifold reproduces N_o = 1112203)
for s in ['water','bmimTf2N','emimTf2N','bmimDCA','bmimSCN']:
    for Ug in [0.01, 0.05, 0.1, 0.5, 1.0, 5.0]:
        rows.append(dict(caseName=f'screen_{s}_{Ug:g}mps', sorbent=s, stage='abs', towerHeight=3, towerDiameter=12,
                         orificeDensity=10000, d_orifice_mm=0.8, Q_orifice=3.75e-7, T_C=20, P=1, P_baseline=1,
                         Ug_baseline=Ug, ppm0_CO2=400, ppm0_O2=210000, ppm0_N2=780000, H_scale=1, Ncol=1))

#BGAL parametric cases
for name in wb_in.sheetnames:
    ws = wb_in[name]
    sorbent = 'bmimDCA' if name.startswith('IL') else 'water'
    tag = name.replace(' - ', '_').replace(' ', '_')
    a = grab(ws, 0, dict(ABS, Ncol='Ncol'))
    d = grab(ws, 4, dict(DES, Ncol='Ncol'))

    a.update(caseName=tag, sorbent=sorbent, stage='abs', H_scale=H_scale.get(name, 1))
    if name in pressureCases:
        a.update(P_baseline=10, Ug_baseline=0.05) #mass flow fixed at the 10 atm / 0.05 m/s baseline
    else:
        a.update(P_baseline=a['P'], Ug_baseline=a['Ug'])
    rows.append(a)

    d.update(caseName=tag, sorbent=sorbent, stage='des', H_scale=H_scale.get(name, 1))
    if name == 'IL0.05desorb': #sheet lists 0.1 atm
        d.update(P=0.05)
    if name == 'IL0.2desorb': #sheet lists 0.1 atm
        d.update(P=0.2)
    rows.append(d)

wb = openpyxl.Workbook()
ws = wb.active
ws.title = 'Cases'
ws.append([n for n, _, _ in cols])
fill = {0: PatternFill('solid', fgColor='E8E8E8'), 2: PatternFill('solid', fgColor='DDEBF7')}
rowOf = {}
for i, r in enumerate(rows, start=2):
    rowOf[(r['caseName'], r['stage'])] = i
    f = formulas(i, r['stage'], r['sorbent'])
    for j, (n, _, kind) in enumerate(cols, start=1):
        cell = ws.cell(row=i, column=j, value=r.get(n) if kind == 1 else f.get(n))
        if kind in fill:
            cell.fill = fill[kind]

#orifice density links (as in the BGAL Settings sheets)
optRow = rowOf[('Water_Opt', 'abs')]
for stage in ['abs', 'des']:
    ws[f'{L["orificeDensity"]}{rowOf[("Water_High_GasLH", stage)]}'] = f'={L["orificeDensity"]}{optRow}*20'

#gas velocity links: High/Low GasLH keep the Opt L/G (gas scaled with liquid, as in the BGAL Settings sheets)
ws[f'{L["Ug_baseline"]}{rowOf[("Water_High_GasLH", "abs")]}'] = f'={L["Ug_baseline"]}{optRow}*20'
ws[f'{L["Ug_baseline"]}{rowOf[("Water_Low_GasLH", "abs")]}'] = f'={L["Ug_baseline"]}{optRow}*0.1'
#IL High/Low Abs are the IL Opt column at a different height scale
ilOptRow = rowOf[('IL_Opt', 'abs')]
for name in ['IL_High_Abs', 'IL_Low_Abs']:
    ws[f'{L["Ug_baseline"]}{rowOf[(name, "abs")]}'] = f'={L["Ug_baseline"]}{ilOptRow}'

for j in range(1, len(cols)+1):
    ws.cell(row=1, column=j).font = Font(bold=True)
ws.freeze_panes = 'D2'

#column units and notes
wc = wb.create_sheet('Columns')
wc.append(['column', 'units', 'type'])
for n, u, kind in cols:
    wc.append([n, u, ['formula', 'input', 'MATLAB output'][kind]])

wb.save(out)
print(len(rows), 'rows ->', out)

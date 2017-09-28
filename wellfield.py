               
import numpy as np
import flopy
import os
import pandas as pd
import math
from get_layerBudget import *
from create_model import calc_model_wellcoordinates


   
print("###############################################################################")
print("Modell Brunnenfeld")
print("###############################################################################\n")

 ### 1) Modell-Erstellung ###

    ### Modell-Domäne

Ly = 5400. #Erstreckung in y-Richung
Lx = 5000. #Erstreckung in x-Richung

#Ausrichtung des Modells in Abhängigkeit der Brunnenkoordinaten
upleft_coord = calc_model_wellcoordinates(Ly = Ly, 
                                          Lx = Lx, 
                                          csvDir = '.',
                                          csvFile = 'wells_nodes.csv',
                                          exp_dir = '.')
xul = upleft_coord['xul']
yul = upleft_coord['yul']
proj4_str = 'ESPG:31466'


    ### Grid-Definition
    
ztop = 200. #Höhe des Modells
zbot = 0. #Basis des Modells
nlay = 3 #Anzahl Schichten
grid_spacing = 50 #Abstand der Gitternetzlinien
delr = grid_spacing #Abstand der Reihen-Gitternetzlinien
delc = grid_spacing #Abstand der Säulen-Gitternetzlinien
delv = np.array([155, 15, 30], dtype=np.float32) #Schichtmächtigkeiten (für jede layer)
nrow = int(Ly / delr) #Anzahl Reihen
ncol = int(Lx / delc) #Anzahl Säulen
start_datetime = '1/1/2007' #Startzeit


 ### 2) Modell-Parametrierung ###
    
#horizontal conductivity
hk = np.array([2e-5*3600*24, 2e-9*3600*24, 3e-5*3600*24],dtype=np.float32) # (für jede layer)
#vertical conductivity
vka =  np.array([2e-5*3600*24, 2e-9*3600*24, 3e-5*3600*24],dtype=np.float32) # (für jede layer)
#specific yield
sy = np.array([0.123, 0.023, 0.15],dtype=np.float32) # (für jede layer)
#specific storage
ss = np.array([1.e-4, 1.e-4, 1.e-4],dtype=np.float32) # (für jede layer)
#layer type
laytyp = np.int_([1, 1, 1]) # 1 - ungespannt, 0 - gespannt (für jede layer)
#Ausgangswasserspiegel
head_north = np.array([160, 160, 115], dtype=np.float32)
#Grundwassergradient
head_gradient_northSouth = 0
#Verkippung der Schichten
gradient_northSouth = 50/Ly #definiere einen linearen Gradient  
#Lekage über Bohrlöcher
kf_borehole = 2e-9*24*3600 #kf-Wert (m/Tag)


    ### Zeitliche Diskretisierung
    
totsim = 12*365 #Gesamtlaufzeit in Tagen
nper = 12 #Anzahl der Zeitschritte


    ### Zielschicht für die Ergebnisdarstellung
    
plot_layer = 2 


    ### Verkippung der Schichten
    
def set_layerbottom(botm_north, 
                    gradient_northSouth
                    ):
    mybotm = np.zeros((nlay, nrow, ncol), dtype=np.float32)
    x = np.linspace(0, 1, ncol)
    for layer, value in enumerate(botm_north):
       elev_south = value+gradient_northSouth*Ly
       print("Setting bottom for layer " + str(layer) + " (elev. north: " + 
             str(botm_north[layer]) + " m, elev. south: " +   str(elev_south) + " m)")
       y = np.linspace(botm_north[layer], elev_south, nrow)
       xv, yv = np.meshgrid(x, y)
       mybotm[layer] = yv
    return(mybotm)
botm = set_layerbottom(botm_north = np.array([ztop - delv[0],ztop - sum(delv[0:2]), zbot], 
                                             dtype=np.float32), 
                       gradient_northSouth = gradient_northSouth #definiere einen linearen Gradient  
                       )

       
       ### Variablen für das BAS-Package
       
#Randbedingungen
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
ibound[2, :, 0] = -1 
ibound[:, :, -1] = 1
ibound[0, :, 0] = -1

#Ausgangswasserspiegel
strt = set_layerbottom(botm_north = head_north, gradient_northSouth = head_gradient_northSouth)

#Einhaltung des von Neumann-Kriteriums
def vonNeumann_max_dt(transmiss ,
                      s,
                      dx,
                      const = 0.49
                      ):

    return(s*dx**2/(2*transmiss/const))
max_dt = min(vonNeumann_max_dt(transmiss = hk*delv, 
                  s = ss, 
                  dx = delr))
print('Setting time step to:', str(max_dt), '(days) for all stress periods')


    ### Zeitliche Diskretisierung
    
steady = np.repeat(False, nper) #type of stress period
t_perPeriod = totsim/nper 
perlen = np.repeat(t_perPeriod,nper) #length of a stress period
nstp =  np.ceil(perlen/1).astype(int) #number of time steps in a stress period


     ### Modflow Packages
    
modelname = 'wellfield'
mf = flopy.modflow.Modflow(modelname)
dis = flopy.modflow.ModflowDis(mf, #model discretisation
                               nlay, 
                               nrow, 
                               ncol, 
                               delr = delr, 
                               delc = delc,
                               top = ztop, 
                               botm = botm,
                               nper = nper, 
                               perlen = perlen, 
                               nstp = nstp,
                               tsmult = 1, 
                               steady = steady,
                               start_datetime = start_datetime)
                               
bas = flopy.modflow.ModflowBas(mf, 
                               ibound = ibound, #boundary conditions
                               strt = strt #starting heads
                               #strt = head[1,:,:]
                               )
lpf = flopy.modflow.ModflowLpf(mf, #layer-property-flow
                               hk = hk, 
                               vka = vka, 
                               sy = sy, 
                               ss = ss, 
                               laytyp = laytyp,
                               constantcv = True)
pcg = flopy.modflow.ModflowPcg(mf,
                               hclose = 5E-4,
                               rclose = 1E-3) #Preconditioned Conjugate-Gradient

                              
    ### Brunnenrandbedingung / MNW-2 package 

node_data  = pd.read_csv('wells_nodes.csv')
node_data["k"] = node_data["k"].astype(int)
node_data["j"] = (node_data["x"]/delc).astype(int)
node_data["i"] = (node_data["y"]/delr).astype(int)
node_data["j"] = (node_data["x"]/delc).astype(int)
wells_location = node_data


    ### Lekage über Bohrlöcher und Brunnen

def get_realLeakage(area_welllocs = 0.3, #meter^2
                    area_model = 2500,  #meter^2 
                    kf_welllocs = 1E-7, #meter/day
                    kf_natural = 1E-6 #meter/day
                    ):
    return((area_welllocs * kf_welllocs + (area_model - area_welllocs) * kf_natural)/area_model);
area_borehole = 0.3 ###meter^2
kf_borehole = kf_borehole #### meter / day
hk_with_boreholes = lpf.hk.array  ###copied from initial model 
vka_with_boreholes = lpf.vka.array ###copied from initial model 

# Replacing natural leakage with combined leakage value for all wells screened 
# below MODFLOW layer 2 (i.e. in flopy: below layer 1) 
for well_cell_idx in np.arange(0,len(wells_location),1):
    tmp_well = wells_location.ix[[well_cell_idx]]
    k = int(tmp_well.ix[:, ['k']].values)
    if (k >= 2):
        i = int(tmp_well.ix[:,['i']].values)
        j = int(tmp_well.ix[:,['j']].values)
        leak_layer = 1
        area_model = dis.delc.array[i]*dis.delr.array[j] - area_borehole
    
        hk_with_boreholes[leak_layer,i,j] = get_realLeakage(area_borehole, area_model, kf_borehole, lpf.hk.array[leak_layer,i,j])
        vka_with_boreholes[leak_layer,i,j] = get_realLeakage(area_borehole, area_model, kf_borehole, lpf.vka.array[leak_layer,i,j])
    else: 
        print('Well not screened in model layer 3 or higher')

#Overwrite existing lpf package with combined natural+borehole leakage for layer 1        
lpf = flopy.modflow.ModflowLpf(mf, #layer-property-flow
                               hk = hk_with_boreholes, 
                               vka = vka_with_boreholes, 
                               sy = sy, 
                               ss = ss, 
                               laytyp = laytyp,
                               constantcv = True)

node_data = node_data.to_records()
node_data


stress_period_data  = pd.read_csv('wells_times.csv')

wells_info =  pd.merge(left=wells_location,
                       right=stress_period_data, 
                       left_on='wellid', 
                        right_on='wellid')
wells_info = wells_info[wells_info['qdes'] != 0]

pers = stress_period_data.groupby('per')
stress_period_data = {i: pers.get_group(i).to_records() for i in range(0,nper)}

nwells = len(node_data)
                      
mnw2 = flopy.modflow.ModflowMnw2(model=mf, mnwmax=nwells,
                 node_data=node_data, 
                 stress_period_data=stress_period_data,
                 itmp=np.repeat(nwells, nper), # reuse second per pumping for last stress period
                 extension='mnw2', 
                 unitnumber=23
                 )


mnw2.write_file(modelname + ".mnw2")


     ### Output control

output_features = ['save head', 
                   'save budget']
each_time_step = True ### if True (for all time steps), if False only for last time 
#for each stress period 
                   
ocdict = {}
for sper in range(0,nper):
    if  steady[sper]==False:
        if each_time_step == True:
            for time_step in np.arange(0, nstp[sper], 1):
                key = (sper, time_step)
                ocdict[key] = output_features
                key = (sper+1, 0)  
                ocdict[key] = []  
        else:
            key = (sper, nstp[sper]-1)
            ocdict[key] = output_features
            key = (sper+1, 0)  
            ocdict[key] = []  
    else:
        key = (sper, 0)
        ocdict[key] = output_features
        key = (sper+1, 0)  
        ocdict[key] = [] 

oc = flopy.modflow.ModflowOc(mf, 
                            stress_period_data = ocdict,
                             compact=True)


     ### Write the model input files

mf.write_input()


    ### Adding MNW2 & MWNI in MODFLOW input file 

if 'mnw2' in locals():
    print("Writing MNW2 & MNWI package to " + modelname + ".nam file")
    with open(modelname + ".nam", 'a') as f:
        f.write('MNW2          23 ' + modelname + '.mnw2')
        f.write('\nMNWI          79 ' + modelname + '.mnwi')
        f.write('\nDATA          82 ' + modelname + '.byn')
else:
    print("Not using MNW2 & MNWI package")

    

 ### 3) Run the model ###

success, mfoutput = mf.run_model(silent=False, pause=False)
if not success:
    raise Exception('MODFLOW did not terminate normally.')

  
    
 ### 4) Import Modell Ergebnisse ###
#
#import matplotlib.pyplot as plt
#import flopy.utils.binaryfile as bf
#
#### Get the headfile and budget file objects
#headobj = bf.HeadFile(modelname+'.hds')
#times = headobj.get_times()
#cbb = bf.CellBudgetFile(modelname+'.cbc')
#
#### Save head file to shape
#headobj.to_shapefile('test_heads_sp12.shp', totim=times[-1], mflay=2, attrib_name='head')


import os
import numpy as np
import flopy
import pandas as pd
import sys
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf
from get_layerBudget import *
from create_model import calc_model_wellcoordinates
from observations import *


 ### 1) Import Modell Ergebnisse ###

### Get the headfile and budget file objects
headobj = bf.HeadFile(modelname+'.hds')
times = headobj.get_times()
cbb = bf.CellBudgetFile(modelname+'.cbc')

### Save head file to shape
headobj.to_shapefile('test_heads_sp12.shp', totim=times[-1], mflay=2, attrib_name='head')

 ### 2) Prozessierung und Darstellung Ergebnisse  ###
    
    ### Wasserbilanz
    
### Get layer based budget for each time step
layer_budget = get_layerbudget("wellfield", 
                               nper, 
                               perlen, 
                               nlay, 
                               debug = True)

### Aggregate budget for layer for whole simulation
layer_budget_perLayer = layer_budget.groupby(['layer']).sum().reset_index()

### Aggregate budget for stress period & layer
layer_budget_perStressPeriod = layer_budget.groupby(['layer', 'stress_period']).sum().reset_index()

### Filter only lowest layer
layer3_budget_perStressPeriod = layer_budget_perStressPeriod[layer_budget_perStressPeriod['layer'] == plot_layer]
mf_list = flopy.utils.MfListBudget(modelname+".list")
budget_incremental, budget_cumulative = mf_list.get_dataframes(start_datetime='31-12-2006')
layer3_budget_perStressPeriod['MNW2_IN'] = np.append(budget_cumulative['MNW2_IN'][0],
                                                    budget_cumulative['MNW2_IN'].diff().as_matrix()[1:])
layer3_budget_perStressPeriod['MNW2_OUT'] = np.append(budget_cumulative['MNW2_OUT'][0],
                                                    budget_cumulative['MNW2_OUT'].diff().as_matrix()[1:])
layer3_budget_perStressPeriod['LEAKAGE_FROM_LAYER2'] = layer_budget_perStressPeriod[layer_budget_perStressPeriod['layer'] == 1]['FLOW_LOWER_FACE'].as_matrix()

### Grafik Wasserbilanz
bar_width = 0.35
#plt.bar(layer3_budget_perStressPeriod['stress_period'], layer3_budget_perStressPeriod['MNW2_OUT'], bar_width, color='b', label="QBrunnen")
#plt.bar(layer3_budget_perStressPeriod['stress_period'], 
#        layer3_budget_perStressPeriod['CONSTANT_HEAD_OUT'], 
#        bar_width, 
#        color= 'green', 
#        label="Out: Rand", 
#        bottom=layer3_budget_perStressPeriod['MNW2_OUT'])
#plt.bar(layer3_budget_perStressPeriod['stress_period'], 
#        layer3_budget_perStressPeriod['STORAGE_OUT'], 
#        bar_width, 
#        color= 'brown', 
#        label="Out: Vorrat", 
#        bottom=layer3_budget_perStressPeriod['MNW2_OUT']+layer3_budget_perStressPeriod['CONSTANT_HEAD_OUT'])
plt.bar(layer3_budget_perStressPeriod['stress_period'] + bar_width,  
        layer3_budget_perStressPeriod['STORAGE_IN'], 
        bar_width, 
        color='r', 
        label="QVorrat")
plt.bar(layer3_budget_perStressPeriod['stress_period'] + bar_width, 
        layer3_budget_perStressPeriod['CONSTANT_HEAD_IN'], 
        bar_width,
        color= 'orange', 
        label="QRand", 
        bottom=layer3_budget_perStressPeriod['STORAGE_IN'])
plt.bar(layer3_budget_perStressPeriod['stress_period'] + bar_width, 
        layer3_budget_perStressPeriod['LEAKAGE_FROM_LAYER2'],
        bar_width,
        color= 'grey', 
        label="QLeakage", 
        bottom=layer3_budget_perStressPeriod['CONSTANT_HEAD_IN'] + layer3_budget_perStressPeriod['STORAGE_IN'])
plt.legend(bbox_to_anchor=(1.2, 0.9), bbox_transform=plt.gcf().transFigure)
plt.title('Entnahme aus 6B (in m3 pro Stressperiode)')
plt.axis([0, 12, 0, 1.8e7])
plt.ylabel('m3')
plt.xlabel('Stressperiode')
plt.xticks(layer3_budget_perStressPeriod['stress_period'])
plt.xticks(layer3_budget_perStressPeriod['stress_period'])
plt.savefig('budget_layer3.png', dpi=300, bbox_inches='tight')
plt.show()


    ### Gleichenpläne für Stressperioden
    
### Setup contour parameters
levels = np.linspace(0, 80, 17)
extent = (delr/2., Ly - delr/2., delc/2., Lx - delc/2.)
print('Levels: ', levels)
print('Extent: ', extent)
timeslist = [t_perPeriod,t_perPeriod*2, t_perPeriod*3, t_perPeriod*4, t_perPeriod*5,t_perPeriod*6,t_perPeriod*7,t_perPeriod*7,t_perPeriod*8,t_perPeriod*9,t_perPeriod*10,t_perPeriod*11,t_perPeriod*12]

### Well point
def getStressPeriod(time, 
                    stress_period_ende = perlen):
    cum_perende = np.cumsum(stress_period_ende)
    stress_period = sum(time > cum_perende)
    return(stress_period )
    
### Make the plots
for iplot, time in enumerate(timeslist):
    print('*****Processing time: ', time)
    
    # Get head data
    head = headobj.get_data(totim=time)
    
    # Get pumping wells 
    str_per = getStressPeriod(time)
 
    # Print statistics
    print('Head statistics')
    print('  min: ', head.min())
    print('  max: ', head.max())
    print('  std: ', head.std())

    # Extract flow right face and flow front face
    frf = cbb.get_data(text='FLOW RIGHT FACE', totim=time)[0]
    fff = cbb.get_data(text='FLOW FRONT FACE', totim=time)[0]

    # Create the plot
    fig = plt.figure(figsize=(8, 8))
    plt.subplot(1, 1, 1, aspect='equal')
    plt.title('Jahr: ' + str(int(time/365+2006)) + ' \n(Stressperiode: ' + str(str_per+1) + ")")
    modelmap = flopy.plot.ModelMap(model=mf, layer=plot_layer)
    qm = modelmap.plot_ibound()
    lc = modelmap.plot_grid()
    cs = modelmap.contour_array(head, levels=levels)
    plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f', zorder=11)
    quiver = modelmap.plot_discharge(frf, fff, head=head)
    mfc='black'
    plt.plot(xkoord_1-xul, Ly-(yul-ykoord_1), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
    plt.plot(xkoord_2-xul, Ly-(yul-ykoord_2), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
    plt.plot(xkoord_3-xul, Ly-(yul-ykoord_3), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
    plt.plot(xkoord_4-xul, Ly-(yul-ykoord_4), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
    plt.plot(xkoord_5-xul, Ly-(yul-ykoord_5), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
    plt.plot(xkoord_6-xul, Ly-(yul-ykoord_6), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
    if 'mnw2' in locals():
        print("Using MNW2 package and plotting active wells")
        mnw_wells = wells_info[wells_info['per'] == str_per] 
        plt.plot(mnw_wells['j']*delc, 
                 Ly-(mnw_wells['i']*delr), 
                 lw=0, 
                 marker='o', 
                 markersize=3, 
                 markeredgewidth=1,
                 markeredgecolor='black', 
                 markerfacecolor=mfc, 
                 zorder=9)
    else:
        print("Using WEL package and plotting active wells")
        pumping_wells = pumping_data[str_per]
        plt.plot(pumping_wells[:,2]*delc, 
                 pumping_wells[:,1]*delr, 
                 lw=0, marker='o', markersize=3, 
                 markeredgewidth=1,
                 markeredgecolor='black', markerfacecolor=mfc, zorder=9)
    plt.savefig('Aktive Brunnen.png', dpi=300, bbox_inches='tight')
    plt.show()

    
    ### Restwassermächtigkeiten und Gleichen zum Endzeitpunkt 
   
head = headobj.get_data(totim=times[len(times)-1])
time = times[len(times)-1]
mytitle = 'GW-Gleichen und Restwassermächtigkeit im 6B in 2015' 
ax = fig.add_subplot(1, 1, 1, aspect='equal')
title = ax.set_title(mytitle)
modelmap = flopy.plot.ModelMap(model=mf, rotation=0)
rw = modelmap.contour_array(head[plot_layer,:,:]-botm[plot_layer], levels=levels, colors= 'brown')
levels = np.linspace(-100, 20, 25)
hh = modelmap.contour_array(head[plot_layer,:,:]-110, levels=levels, colors= 'blue', linestyles='solid')
plt.clabel(rw, inline=1, fontsize=10, fmt='%1.f', zorder=11)
plt.clabel(hh, inline=1, fontsize=10, fmt='%1.f', zorder=11)
plt.plot(xkoord_1-xul, Ly-(yul-ykoord_1), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(xkoord_2-xul, Ly-(yul-ykoord_2), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(xkoord_50494-xul, Ly-(yul-ykoord_50494), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(xkoord_3-xul, Ly-(yul-ykoord_3), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(xkoord_4-xul, Ly-(yul-ykoord_4), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(xkoord_5-xul, Ly-(yul-ykoord_5), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(xkoord_6-xul, Ly-(yul-ykoord_6), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
if 'mnw2' in locals():
        print("Using MNW2 package and plotting active wells")
        mnw_wells = wells_info[wells_info['per'] == sper] 
        plt.plot(mnw_wells['j']*delc, 
                 Ly-(mnw_wells['i']*delr), 
                 lw=0, 
                 marker='o', 
                 markersize=3, 
                 markeredgewidth=1,
                 markeredgecolor='black', 
                 markerfacecolor='black', 
                 zorder=9)
plt.savefig('Gleichen_Restwasser_6B.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()


    ### Grundwasserstands-Zeitreihen

### Lage der Messstellen im Modell    
obsPoint_1_6B = [plot_layer, xkoord_1-xul, yul-ykoord_1]
obsPoint_1_6D = [0, xkoord_1-xul, yul-ykoord_1]
obsPoint_2_6B = [plot_layer, xkoord_2-xul, yul-ykoord_2]
obsPoint_2_6D = [0, xkoord_2-xul, yul-ykoord_2]
obsPoint_3_6B = [plot_layer, xkoord_3-xul, yul-ykoord_3]
obsPoint_3_6D = [0, xkoord_3-xul, yul-ykoord_3]
obsPoint_4_6B = [plot_layer, xkoord_4-xul, yul-ykoord_4]
obsPoint_4_6D = [0, xkoord_4-xul, yul-ykoord_4]
obsPoint_5_6B = [plot_layer, xkoord_5-xul, yul-ykoord_5]
obsPoint_5_6D = [0, xkoord_5-xul, yul-ykoord_5]
obsPoint_6_6B = [plot_layer, xkoord_6-xul, yul-ykoord_6]
obsPoint_6_6D = [0, xkoord_6-xul, yul-ykoord_6]

### Import Messpunkt 1
### Convert observation point to layer, column, row format
idx1 = (obsPoint_1_6B[0], 
       round(obsPoint_1_6B[2]/delr,0), 
       round(obsPoint_1_6B[1]/delc,0))
idx1_1 = (obsPoint_1_6D[0], 
       round(obsPoint_1_6D[2]/delr,0), 
       round(obsPoint_1_6D[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx1[1]), int(idx1[2])]

### Grafik
ts = headobj.get_ts(idx1) 
ts_1 = headobj.get_ts(idx1_1) 
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(obsPoint_1_6B[0] + 1, int(obsPoint_1_6B[1]), int(obsPoint_1_6B[2]))
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_1_6B[:, 0], obs_measured_1_6B[:, 1] + botm[nlay-1,int(idx1[1]),int(idx1[2])], color="blue", ls=':', label='Pegel 6B (814193)')
plt.plot(obs_measured_1_6D[:, 0], obs_measured_1_6D[:, 1] + botm[nlay-1,int(idx1_1[1]),int(idx1_1[2])], color="red", ls=':', label='Pegel 6D (814192)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('81419.png', dpi=300, bbox_inches='tight')
plt.show()


### Import Messpunkt 2
### Convert observation point to layer, column, row format
idx2 = (obsPoint_2_6B[0], 
       round(obsPoint_2_6B[2]/delr,0), 
       round(obsPoint_2_6B[1]/delc,0))
idx2_1 = (obsPoint_2_6D[0], 
       round(obsPoint_2_6D[2]/delr,0), 
       round(obsPoint_2_6D[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx2[1]), int(idx2[2])]

### Grafik                         
ts = headobj.get_ts(idx2)
ts_1 = headobj.get_ts(idx2_1)  
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(obsPoint_2_6B[0] + 1, int(obsPoint_2_6B[1]), int(obsPoint_2_6B[2]))
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_2_6B[:, 0], obs_measured_2_6B[:, 1] + botm[nlay-1,int(idx2[1]),int(idx2[2])], ls=':', label='Pegel 6B (503372)')
plt.plot(obs_measured_2_6D[:, 0], obs_measured_2_6D[:, 1] + botm[nlay-1,int(idx2_1[1]),int(idx2_1[2])], color="red", ls=':', label='Pegel 6D (503371)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('50337.png', dpi=300, bbox_inches='tight')
plt.show()

### Import Messpunkt 3
### Convert observation point to layer, column, row format
idx3 = (obsPoint_3_6B[0], 
       round(obsPoint_3_6B[2]/delr,0), 
       round(obsPoint_3_6B[1]/delc,0))
idx3_1 = (obsPoint_3_6D[0], 
       round(obsPoint_3_6D[2]/delr,0), 
       round(obsPoint_3_6D[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx2[1]), int(idx2[2])]

### Grafik
ts = headobj.get_ts(idx3)
ts_1 = headobj.get_ts(idx3_1)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(obsPoint_3_6B[0] + 1, int(obsPoint_3_6B[1]), int(obsPoint_3_6B[2]))
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_3_6B[:, 0], obs_measured_3_6B[:, 1] + botm[nlay-1,int(idx3[1]),int(idx3[2])], color='b', ls=':', label='Pegel 6B (504943)')
#plt.plot(obs_measured_3_6D[:, 0], obs_measured_3_6D[:, 1] + botm[nlay-1,int(idx3_1[1]),int(idx3_1[2])], color='r', ls=':', label='Pegel 6D (504943)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('50261.png', dpi=300, bbox_inches='tight')
plt.show()


### Import Messpunkt 4
obs_measured_503762  = np.loadtxt('obs_head_time_503762.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_503761  = np.loadtxt('obs_head_time_503761.csv', 
                           delimiter=",",
                           skiprows = 1)

### User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_503762 = [plot_layer, xkoord_50376-xul, yul-ykoord_50376]
obsPoint_503761 = [0, xkoord_50376-xul, yul-ykoord_50376]

### Convert observation point to layer, column, row format
idx5 = (obsPoint_503762[0], 
       round(obsPoint_503762[2]/delr,0), 
       round(obsPoint_503762[1]/delc,0))
idx5_1 = (obsPoint_503761[0], 
       round(obsPoint_503761[2]/delr,0), 
       round(obsPoint_503761[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx5[1]), int(idx5[2])]

### Grafik                                                  
ts = headobj.get_ts(idx5)
ts_1 = headobj.get_ts(idx5_1)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(obsPoint_503762[0] + 1, int(obsPoint_503762[1]), int(obsPoint_503762[2]))
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_503762[:, 0], obs_measured_503762[:, 1] + botm[nlay-1,int(idx5[1]),int(idx5[2])], ls=':', label='Pegel 6B (503762)')
plt.plot(obs_measured_503761[:, 0], obs_measured_503761[:, 1] + botm[nlay-1,int(idx5_1[1]),int(idx5_1[2])], color="red", ls=':', label='Pegel 6D (503761)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('50376.png', dpi=300, bbox_inches='tight')
plt.show()

### Import Messpunkt 5
obs_measured_805722  = np.loadtxt('obs_head_time_805722.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_805721  = np.loadtxt('obs_head_time_805721.csv', 
                           delimiter=",",
                           skiprows = 1)

### User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_805722 = [plot_layer, xkoord_80572-xul, yul-ykoord_80572]
obsPoint_805721 = [0, xkoord_80572-xul, yul-ykoord_80572]

### Convert observation point to layer, column, row format
idx6 = (obsPoint_805722[0], 
       round(obsPoint_805722[2]/delr,0), 
       round(obsPoint_805722[1]/delc,0))
idx6_1 = (obsPoint_805721[0], 
       round(obsPoint_805721[2]/delr,0), 
       round(obsPoint_805721[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx6[1]), int(idx6[2])]

### Grafik                                        
ts = headobj.get_ts(idx6)
ts_1 = headobj.get_ts(idx6_1)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(obsPoint_805722[0] + 1, int(obsPoint_805722[1]), int(Ly-obsPoint_805722[2]))
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_805722[:, 0], obs_measured_805722[:, 1] + botm[nlay-1,int(idx6[1]),int(idx6[2])], ls=':', label='Pegel 6B (805722)')
plt.plot(obs_measured_805721[:, 0], obs_measured_805721[:, 1] + botm[nlay-1,int(idx6_1[1]),int(idx6_1[2])], color="red", ls=':', label='Pegel 6D (805721)')
plt.axhspan(ymin=idx_bot[1], ymax=idx_bot[0], color='lightgrey', label='6C')
plt.axhspan(ymin=idx_bot[2], ymax=idx_bot[1], color='bisque', label='6B')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, idx_bot[2], idx_bot[0]])
plt.savefig('80572.png', dpi=300, bbox_inches='tight')
plt.show()


### Import Ziellinie 1
obs_target_504943  = np.loadtxt('obs_head_time_504943.csv', 
                           delimiter=",",
                           skiprows = 1)

### User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_504943 = [plot_layer, xkoord_50494-xul, yul-ykoord_50494]

### Convert observation point to layer, column, row format
idx = (obsPoint_504943[0], 
       round(obsPoint_504943[2]/delr,0), 
       round(obsPoint_504943[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx[1]), int(idx[2])]
                         
### Grafik                                        
ts = headobj.get_ts(idx)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(int(obsPoint_504943[0]) + 1, int(obsPoint_504943[1]), int(Ly-obsPoint_504943[2]))
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(obs_target_504943[:, 0], obs_target_504943[:, 1] + botm[nlay-1,int(idx[1]),int(idx[2])], ls=':', label='Ziellinie 6B (504943)')
plt.axhspan(ymin=idx_bot[1], ymax=idx_bot[0], color='lightgrey', label='6C')
plt.axhspan(ymin=idx_bot[2], ymax=idx_bot[1], color='bisque', label='6B')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, idx_bot[2], idx_bot[0]])
plt.savefig('Zielinie_504943', dpi=300, bbox_inches='tight')
plt.show()

### Import Ziellinie 2
obs_target_504952  = np.loadtxt('target_head_time_504952.csv', 
                           delimiter=",",
                           skiprows = 1)
### Import Pegel
obs_measured_504952  = np.loadtxt('obs_head_time_504952.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_504952 = [plot_layer, xkoord_50495-xul, yul-ykoord_50495]
## Convert observation point to layer, column, row format
idx = (obsPoint_504952[0], 
       round(obsPoint_504952[2]/delr,0), 
       round(obsPoint_504952[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx[1]), int(idx[2])]
                         
### Grafik                         
ts = headobj.get_ts(idx)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(int(obsPoint_504952[0]) + 1, int(obsPoint_504952[1]), int(Ly-obsPoint_504952[2]))
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(obs_target_504952[:, 0], obs_target_504952[:, 1] + botm[nlay-1,int(idx[1]),int(idx[2])], ls=':', color="blue", label='Ziellinie 6B (504952)')
plt.plot(obs_measured_504952[:, 0], obs_measured_504952[:, 1] + botm[nlay-1,int(idx[1]),int(idx[2])], ls=':', color="red",label='Pegel 6B (504952)')
plt.axhspan(ymin=idx_bot[1], ymax=idx_bot[0], color='lightgrey', label='6C')
plt.axhspan(ymin=idx_bot[2], ymax=idx_bot[1], color='bisque', label='6B')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 5500, idx_bot[2], idx_bot[0]])
plt.savefig('Zielinie_504952.png', dpi=300, bbox_inches='tight',transparent=True)
plt.show()


### Modell Profilschnitte

### West-Ost-Schnitt
modelxsect = flopy.plot.ModelCrossSection(model=mf, line={'row': 80}) #Eingabe Modellreihe
patches = modelxsect.plot_ibound(head=head)
wt = modelxsect.plot_surface(head, masked_values=[999.], color='blue', lw=1)
linecollection = modelxsect.plot_grid()
plt.title('Profilschnitt in W-O-Richtung mit Grundwasserständen')
plt.ylabel('m')
plt.xlabel('m')
plt.axis([0, 5000, 0, 180])
plt.savefig('xsect.png', dpi=300, bbox_inches='tight')
plt.show()


### Nord-Süd-Schnitt
modelysect = flopy.plot.ModelCrossSection(model=mf, line={'column': 70}) #Eingabe Modellspalte
patches = modelysect.plot_ibound(head=head)
wt = modelysect.plot_surface(head, masked_values=[999.], color='blue', lw=1)
linecollection = modelysect.plot_grid()
plt.title('Profilschnitt in N-S-Richtung mit Grundwasserständen')
plt.ylabel('m')
plt.xlabel('m')
plt.axis([0, 5400, 0, 200])
plt.savefig('ysect.png', dpi=300, bbox_inches='tight')
plt.show()

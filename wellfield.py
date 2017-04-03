import numpy as np
import flopy
import os
import pandas as pd
import math
from get_layerBudget import *

# Model domain and grid definition
Ly = 5400.
Lx = 5000.
ztop = 200.
zbot = 0.
nlay = 3
grid_spacing = 50
delr = grid_spacing
delc = grid_spacing
delv = np.array([160, 15, 25], dtype=np.float32)
nrow = int(Ly / delr)
ncol = int(Lx / delc)
start_datetime = '1/1/2007'
xoff = 1000

###layer decline
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
                       gradient_northSouth = 65/Ly #0.012037037037037037  
                       )

#botm = np.array([ztop - delv[0],ztop - sum(delv[0:2]), zbot], dtype=np.float32)


hk = np.array([2e-5*3600*24, 5e-10*3600*24, 5e-5*3600*24], #horizontal conductivity
              dtype=np.float32)
vka =  np.array([2e-5*3600*24, 5e-10*3600*24, 5e-5*3600*24], #vertical conductivity
                dtype=np.float32)
sy = np.array([0.123, 0.023, 0.15], #specific yield
              dtype=np.float32)
ss = np.array([1.e-4, 1.e-4, 1.e-4], #specific storage
              dtype=np.float32)
laytyp = np.int_([1, 1, 1]) # 1 - ungespannt, 0 - gespannt

#plot_layer = nlay-1

# Variables for the BAS package
# Note that changes from the previous tutorial!
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
ibound[2, :, 0] = -1
ibound[:, :, -1] = 1
ibound[0, :, 0] = -1
#ibound[0:2, 0, 0:15] = -1
#ibound[0:2, -1, 0:15] = -1
#starting heads
head_north = np.array([160, 160, 115], dtype=np.float32)
head_gradient_northSouth = -35/5400

#iniHead= 110
#strt = iniHead * np.ones((nlay, nrow, ncol), dtype=np.float32)
#strt[0, :, :] = 120
strt = set_layerbottom(botm_north = head_north, gradient_northSouth = head_gradient_northSouth)

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

### Desired total simulation time
totsim = 12*365


# Time step parameters
nper = 12 #number of stress periods
steady = np.repeat(False, nper) #type of sress period
t_perPeriod = totsim/nper 
perlen = np.repeat(t_perPeriod,nper) #length of a stress period
#nstp = np.ceil(perlen[1:nper]/max_dt).astype(int)
nstp =  np.ceil(perlen/1).astype(int) #number of time steps in a stress period

# Flopy objects
modelname = 'wellfield'
#mf = flopy.modflow.Modflow(modelname, 
#                           version = 'mf2k',
#                           exe_name='mf2k')
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

# Make list for stress period 1
stageleft = 70. #head on the left boundary
# stageright = 80  #head on the right boundary
bound_sp1 = [] #boundary conditions for stress period 1
for il in range(nlay):
    condleft = hk[il] * delv[il] #conductance on the left boundary
#    condright = hk[il] * (stageright - zbot) * delc #conductance on the right boundary
    for ir in range(nrow):
        bound_sp1.append([il, ir, 0, stageleft, condleft])
#        bound_sp1.append([il, ir, ncol - 1, stageright, condright])
#print('Adding ', len(bound_sp1), 'GHBs for stress period 1.')

#boundary_data = {0: bound_sp1}

# Create the flopy ghb object
#ghb = flopy.modflow.ModflowGhb(mf, stress_period_data = boundary_data)

# Create the well package
# Remember to use zero-based layer, row, column indices!
#
#def create_wellfield(layers = [2],
#                     spacing_x = 1000, 
#                     spacing_y = 1000, 
#                     extent_x = [2600,4800],
#                     extent_y = [500,4500],
#                     pumping_rate = 0):
#
#    layers = np.array(layers)
#    extent_x = np.array(extent_x) #m
#    extent_y = np.array(extent_y) #m
#
#    rows = np.array(range(extent_y[0],
#                      extent_y[1]+1,
#                      spacing_y)) / delr
#                       
#    rows = np.round(rows,0)
#
#    cols = np.array(range(extent_x[0],
#                       extent_x[1]+1,
#                       spacing_x)) / delc
#                       
#    cols = np.round(cols,0)
#    num_wells = rows.size * cols.size
#    wellfield = pd.DataFrame()
#
#    for layer in layers:
#        for well_row in rows:
#            for well_col in cols:
#                tmp = pd.DataFrame([[layer,well_row, well_col, pumping_rate]],
#                               columns=['layer','row','column','pumping_rate'])
#                wellfield = wellfield.append(tmp,ignore_index=True)
#    msg = 'Added ' +  str(num_wells) + ' pumping wells in layer ' + str(layer) + ' each with Q = ' + str(pumping_rate) + ' m3/day' 
#    print(msg)
#    return(wellfield)
#wellfield_df = pd.DataFrame(wellfield,columns=["layer","row","column","pumping_rate"])              
#pd.DataFrame.as_matrix(wellfield_df)                  


## Well package 
#wel_sp1 = create_wellfield()
#wel_sp3 = create_wellfield(spacing_x = 200, spacing_y = 200, pumping_rate = -50)
#wel_sp4 = create_wellfield(spacing_x = 200, spacing_y = 200, pumping_rate = -75)
#wel_sp5 = create_wellfield(spacing_x = 200, spacing_y = 200, pumping_rate = -100)
#wel_sp6 = create_wellfield(spacing_x =200, spacing_y = 200, pumping_rate = -150)
#wel_sp7 = create_wellfield(spacing_x =200, spacing_y = 200, pumping_rate = -200)
#wel_sp8 = create_wellfield(spacing_x =200, spacing_y = 200, pumping_rate = -150)
#wel_sp9 = create_wellfield(spacing_x =200, spacing_y = 200, pumping_rate = -300)                                           
#wel_sp10 = create_wellfield(spacing_x =200, spacing_y = 200,pumping_rate= -800)                                            
#wel_sp11 = create_wellfield(spacing_x =150, spacing_y = 150,pumping_rate=-150)
#wel_sp12 = create_wellfield(spacing_x =150, spacing_y = 150,pumping_rate=-70)
#                                            
#wfs = wel_sp6.append(create_wellfield(spacing_x =1000, spacing_y = 1000, 
#                        extent_x = [2500,2600],
#                        extent_y = [500,4500],
#                        pumping_rate=-50*24))
#
#
#pumping_data =       {0: pd.DataFrame.as_matrix(wel_sp1), 
#                      1: pd.DataFrame.as_matrix(wel_sp3), 
#                      2: pd.DataFrame.as_matrix(wel_sp4),
#                      3: pd.DataFrame.as_matrix(wel_sp5),
#                      4: pd.DataFrame.as_matrix(wel_sp6),
#                      5: pd.DataFrame.as_matrix(wel_sp7), 
#                      6: pd.DataFrame.as_matrix(wel_sp7),
#                      7: pd.DataFrame.as_matrix(wel_sp8),
#                      8: pd.DataFrame.as_matrix(wel_sp8),
#                      9: pd.DataFrame.as_matrix(wel_sp8)}
#                      
#wel = flopy.modflow.ModflowWel(mf, stress_period_data = pumping_data)

### MNW-2 package

node_data  = pd.read_csv('wells_nodes.csv')

node_data["k"] = node_data["k"].astype(int)
node_data["j"] = (node_data["x"]/delc).astype(int)
node_data["i"] = (node_data["y"]/delr).astype(int)
node_data["j"] = (node_data["x"]/delc).astype(int)

wells_location = node_data



### Taking into account higher leakage through boreholes
def get_realLeakage(area_welllocs = 0.3, #meter^2
                    area_model = 2500,  #meter^2 
                    kf_welllocs = 1E-7, #meter/day
                    kf_natural = 1E-6 #meter/day
                    ):
    return((area_welllocs * kf_welllocs + (area_model - area_welllocs) * kf_natural)/area_model);


area_borehole = 0.3 ###meter^2
kf_borehole = 2e-4*24*3600 #### meter / day
hk_with_boreholes = lpf.hk.array  ###copied from initial model 
vka_with_boreholes = lpf.vka.array ###copied from initial model 

### Replacing natural leakage with combined leakage value for all wells screened 
### below MODFLOW layer 2 (i.e. in flopy: below layer 1) 
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
    

        
###Overwrite existing lpf package with combined natural+borehole leakage for layer 1        
lpf = flopy.modflow.ModflowLpf(mf, #layer-property-flow
                               hk = hk_with_boreholes, 
                               vka = vka_with_boreholes, 
                               sy = sy, 
                               ss = ss, 
                               laytyp = laytyp,
                               constantcv = True)
#ids = np.arange(80, 80+node_data["wellid"].count()).astype(str)
#np.array(map(str, ids))
#"DATA          " + ids + " " + node_data["wellid"].values.astype(str) + ".byn"

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
#mf.add_package(mnw2)

#
#
#mnwi = flopy.modflow.mfmnwi.ModflowMnwi(model=mf, 
#                                        wel1flag=0, 
#                                        qsumflag=0, 
#                                        byndflag=82,
#                                        mnwobs = 2,
#                                        wellid_unit_qndflag_qhbflag_concflag = [['well1', 45, 0, 0],
#                                                                                ['well2', 45, 0, 0]],
#                                        unitnumber=79,
#                                        extension='mnwi')



# Output control

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

# Write the model input files
mf.write_input()

###Adding MNW2 & MWNI in MODFLOW input file 
if 'mnw2' in locals():
    print("Writing MNW2 & MNWI package to " + modelname + ".nam file")
    with open(modelname + ".nam", 'a') as f:
        f.write('MNW2          23 ' + modelname + '.mnw2')
        f.write('\nMNWI          79 ' + modelname + '.mnwi')
        f.write('\nDATA          82 ' + modelname + '.byn')
else:
    print("Not using MNW2 & MNWI package")



## Export model data as shapefile
#mf.lpf.hk.export(os.path.join('hk.shp'))
#mf.export(os.path.join('model.shp'))

# Run the model
success, mfoutput = mf.run_model(silent=False, pause=False)
if not success:
    raise Exception('MODFLOW did not terminate normally.')

    
plot_layer = 2 

# Imports
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

# Get the headfile and budget file objects
headobj = bf.HeadFile(modelname+'.hds')
times = headobj.get_times()
cbb = bf.CellBudgetFile(modelname+'.cbc')

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
layer3_budget_perStressPeriod = layer_budget_perStressPeriod[layer_budget_perStressPeriod['layer'] == 2]

mf_list = flopy.utils.MfListBudget(modelname+".list")

budget_incremental, budget_cumulative = mf_list.get_dataframes(start_datetime='31-12-2006')

layer3_budget_perStressPeriod['MNW2_IN'] = np.append(budget_cumulative['MNW2_IN'][0],
                                                    budget_cumulative['MNW2_IN'].diff().as_matrix()[1:])

layer3_budget_perStressPeriod['MNW2_OUT'] = np.append(budget_cumulative['MNW2_OUT'][0],
                                                    budget_cumulative['MNW2_OUT'].diff().as_matrix()[1:])


layer3_budget_perStressPeriod['LEAKAGE_FROM_LAYER2'] = layer_budget_perStressPeriod[layer_budget_perStressPeriod['layer'] == 1]['FLOW_LOWER_FACE'].as_matrix()
 

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
plt.axis([0, 9, 0, 1.8e7])
plt.ylabel('m3')
plt.xlabel('Stressperiode')
plt.xticks(layer3_budget_perStressPeriod['stress_period'])
plt.xticks(layer3_budget_perStressPeriod['stress_period'])
plt.savefig('budget_layer3.png', dpi=300, bbox_inches='tight')
plt.show()


# Setup contour parameters
levels = np.linspace(0, 80, 17)
extent = (delr/2., Ly - delr/2., delc/2., Lx - delc/2.)
print('Levels: ', levels)
print('Extent: ', extent)

# Well point

def getStressPeriod(time, 
                    stress_period_ende = perlen):
    
    cum_perende = np.cumsum(stress_period_ende)
    
    stress_period = sum(time > cum_perende)
   
    return(stress_period )
    

# Make the plots
for iplot, time in enumerate(times):
    print('*****Processing time: ', time)
    
    ### Get head data
    head = headobj.get_data(totim=time)
    
    ### Get pumping wells 
    str_per = getStressPeriod(time)
 
    #Print statistics
    print('Head statistics')
    print('  min: ', head.min())
    print('  max: ', head.max())
    print('  std: ', head.std())

    # Extract flow right face and flow front face
    frf = cbb.get_data(text='FLOW RIGHT FACE', totim=time)[0]
    fff = cbb.get_data(text='FLOW FRONT FACE', totim=time)[0]

    #Create the plot
    #plt.subplot(1, len(mytimes), iplot + 1, aspect='equal')
    plt.subplot(1, 1, 1, aspect='equal')
    plt.title('total time: ' + str(time) + ' days\n(Stress period: ' + str(str_per) + ")")
    modelmap = flopy.plot.ModelMap(model=mf, layer=plot_layer)
    qm = modelmap.plot_ibound()
    lc = modelmap.plot_grid()
    #qm = modelmap.plot_bc('GHB', alpha=0.5)
    cs = modelmap.contour_array(head, levels=levels)
    plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f', zorder=11)
    quiver = modelmap.plot_discharge(frf, fff, head=head)
#    mfc = 'None'
#    if (iplot+1) == len(times):
    mfc='black'
    plt.plot(2050,4950, 
                 lw=0, 
                 marker='o', 
                 markersize=3, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor=mfc, 
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
    plt.show()
    #plt.text(wpt[0]+25, wpt[1]-25, 'well', size=12, zorder=12)


plt.show()



head = headobj.get_data(totim=times[len(times)-1])

time = times[len(times)-1]
#mytitle = 'Heads in layer ' + str(plot_layer) + ' after '+ str(time) + ' days of simulation'
#fig = plt.figure(figsize=(10, 10))
#ax = fig.add_subplot(1, 1, 1, aspect='equal')
#title = ax.set_title(mytitle)
#modelmap = flopy.plot.ModelMap(model=mf, rotation=0)
#quadmesh = modelmap.plot_ibound()
#contour_set = modelmap.plot_array(head[plot_layer,:,:], 
#                                  masked_values=[-1e+30], 
#                                  alpha=0.5)
#cs = modelmap.contour_array(head, levels=levels)
#plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f', zorder=11)
#linecollection = modelmap.plot_grid()
##cb = plt.colorbar(contour_set, shrink=0.4)
#mfc = 'None'
#plt.plot(2050,4950, 
#                 lw=0, 
#                 marker='o', 
#                 markersize=3, 
#                 markeredgewidth=1,
#                 markeredgecolor='red', 
#                 markerfacecolor=mfc, 
#                 zorder=9)
#plt.plot(2650,3250, 
#                 lw=0, 
#                 marker='o', 
#                 markersize=3, 
#                 markeredgewidth=1,
#                 markeredgecolor='red', 
#                 markerfacecolor=mfc, 
#                 zorder=9)
#plt.plot(3900,270, 
#                 lw=0, 
#                 marker='o', 
#                 markersize=3, 
#                 markeredgewidth=1,
#                 markeredgecolor='red', 
#                 markerfacecolor=mfc, 
#                 zorder=9)
#plt.plot(2875,5000, 
#                 lw=0, 
#                 marker='o', 
#                 markersize=3, 
#                 markeredgewidth=1,
#                 markeredgecolor='red', 
#                 markerfacecolor=mfc, 
#                 zorder=9)
#if 'mnw2' in locals():
#        print("Using MNW2 package and plotting active wells")
#        mnw_wells = wells_info[wells_info['per'] == sper] 
#        plt.plot(mnw_wells['j']*delc, 
#                 Ly-(mnw_wells['i']*delr), 
#                 lw=0, 
#                 marker='o', 
#                 markersize=3, 
#                 markeredgewidth=1,
#                 markeredgecolor='black', 
#                 markerfacecolor=mfc, 
#                 zorder=9)
#plt.show()

xul = 2528190
yul = 5662180

##Restwasser 6B
mytitle = 'GW-Gleichen und Restwassermächtigkeit im 6B in 2015' # + str(plot_layer) + ' bottom elevation after '+ str(time) + ' days of simulation'
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
title = ax.set_title(mytitle)
modelmap = flopy.plot.ModelMap(model=mf, rotation=0)
#quadmesh = modelmap.plot_ibound()
#contour_set = modelmap.plot_array(head[plot_layer,:,:]-botm[plot_layer], 
#                                  masked_values=[-1e+30], 
#                                  alpha=0.5)
#levels = np.linspace(0, 80, 17)
#extent = (delr/2., Ly - delr/2., delc/2., Lx - delc/2.)
#print('Levels: ', levels)
#print('Extent: ', extent)
rw = modelmap.contour_array(head[plot_layer,:,:]-botm[plot_layer], levels=levels, colors= 'brown')
levels = np.linspace(-100, 20, 25)
hh = modelmap.contour_array(head[plot_layer,:,:]-110, levels=levels, colors= 'blue', linestyles='solid')
plt.clabel(rw, inline=1, fontsize=10, fmt='%1.f', zorder=11)
plt.clabel(hh, inline=1, fontsize=10, fmt='%1.f', zorder=11)
#linecollection = modelmap.plot_grid()
#cb = plt.colorbar(contour_set, shrink=0.4)
plt.plot(2050+xoff,4950, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(2650+xoff,3250, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(3900+xoff,270, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(2875+xoff,5000, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(1400+xoff,3500, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(1900+xoff,2000, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(2529636-xul+xoff,Ly-(yul-5658103), 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(2530385-xul+xoff,Ly-(yul-5658172), 
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
plt.savefig('Gleichen_Restwasser_6B.png', dpi=300, bbox_inches='tight')
plt.show()

##Restwasser 6D
mytitle = 'GW-Gleichen und Restwassermächtigkeit im 6D in 2015' # + str(plot_layer) + ' bottom elevation after '+ str(time) + ' days of simulation'
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
title = ax.set_title(mytitle)
modelmap = flopy.plot.ModelMap(model=mf, rotation=0)
#quadmesh = modelmap.plot_ibound()
#contour_set = modelmap.plot_array(head[plot_layer,:,:]-botm[plot_layer], 
#                                  masked_values=[-1e+30], 
#                                  alpha=0.5)
#levels = np.linspace(0, 80, 17)
#extent = (delr/2., Ly - delr/2., delc/2., Lx - delc/2.)
#print('Levels: ', levels)
#print('Extent: ', extent)
rw = modelmap.contour_array(head[0,:,:]-botm[0], levels=levels, colors= 'brown')
levels = np.linspace(-70, 50, 25)
hh = modelmap.contour_array(head[0,:,:]-110, levels=levels, colors= 'blue', linestyles='solid')
plt.clabel(rw, inline=1, fontsize=10, fmt='%1.f', zorder=11)
plt.clabel(hh, inline=1, fontsize=10, fmt='%1.f', zorder=11)
#linecollection = modelmap.plot_grid()
#cb = plt.colorbar(contour_set, shrink=0.4)
plt.plot(2050+xoff,4950, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(2650+xoff,3250, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(3900+xoff,270, 
                 lw=0, 
                 marker='o', 
                 markersize=5, 
                 markeredgewidth=1,
                 markeredgecolor='red', 
                 markerfacecolor='red', 
                 zorder=9)
plt.plot(2875+xoff,5000, 
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
plt.savefig('Gleichen_Restwasser_6D.png', dpi=300, bbox_inches='tight')
plt.show()
###Plot the head versus time

### Import measured observation point 1
obs_measured_814193  = np.loadtxt('obs_head_time_814193.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_814192  = np.loadtxt('obs_head_time_814192.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_814193 = [plot_layer, 2050+xoff, 450]
obsPoint_814192 = [0, 2050+xoff, 450]

## Convert observation point to layer, column, row format
idx1 = (obsPoint_814193[0], 
       round(obsPoint_814193[2]/delr,0), 
       round(obsPoint_814193[1]/delc,0))
idx1_1 = (obsPoint_814192[0], 
       round(obsPoint_814192[2]/delr,0), 
       round(obsPoint_814192[1]/delc,0))

idx_bot = dis.botm.array[:, int(idx1[1]), int(idx1[2])]

ts = headobj.get_ts(idx1) 
ts_1 = headobj.get_ts(idx1_1) 
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(obsPoint_814193[0] + 1, obsPoint_814193[1], obsPoint_814193[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_814193[:, 0], obs_measured_814193[:, 1] + botm[nlay-1,int(idx1[1]),int(idx1[2])], color="blue", ls=':', label='Pegel 6B (814193)')
plt.plot(obs_measured_814192[:, 0], obs_measured_814192[:, 1] + botm[nlay-1,int(idx1_1[1]),int(idx1_1[2])], color="red", ls=':', label='Pegel 6D (814192)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('time_series_north.png', dpi=300, bbox_inches='tight')
plt.show()

### Import measured observation point 2
obs_measured_503372  = np.loadtxt('obs_head_time_503372.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_503371  = np.loadtxt('obs_head_time_503371.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_503372 = [plot_layer, 2875+xoff, 400]
obsPoint_503371 = [0, 2875+xoff, 400]
## Convert observation point to layer, column, row format
idx2 = (obsPoint_503372[0], 
       round(obsPoint_503372[2]/delr,0), 
       round(obsPoint_503372[1]/delc,0))
idx2_1 = (obsPoint_503371[0], 
       round(obsPoint_503371[2]/delr,0), 
       round(obsPoint_503371[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx2[1]), int(idx2[2])]
                         
ts = headobj.get_ts(idx2)
ts_1 = headobj.get_ts(idx2_1)  
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = 2875 m and y = 400 m'.format(obsPoint_503372[0] + 1, obsPoint_503372[1], obsPoint_503372[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_503372[:, 0], obs_measured_503372[:, 1] + botm[nlay-1,int(idx2[1]),int(idx2[2])], ls=':', label='Pegel 6B (503372)')
plt.plot(obs_measured_503371[:, 0], obs_measured_503371[:, 1] + botm[nlay-1,int(idx1_1[1]),int(idx1_1[2])], color="red", ls=':', label='Pegel 6D (503371)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('time_series_east.png', dpi=300, bbox_inches='tight')
plt.show()

### Import measured observation point 3
obs_measured_502612  = np.loadtxt('obs_head_time_502612.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_502611  = np.loadtxt('obs_head_time_502611.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_502612 = [plot_layer, 2650+xoff, 2230]
obsPoint_502611 = [0, 2650+xoff, 2230]
## Convert observation point to layer, column, row format
idx2 = (obsPoint_502612[0], 
       round(obsPoint_502612[2]/delr,0), 
       round(obsPoint_502612[1]/delc,0))
idx2_2 = (obsPoint_502611[0], 
       round(obsPoint_502611[2]/delr,0), 
       round(obsPoint_502611[1]/delc,0))
idx_bot = dis.botm.array[:, int(idx2[1]), int(idx2[2])]
ts = headobj.get_ts(idx2)
ts_1 = headobj.get_ts(idx2_2)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = 2650 m and y = 2230 m'.format(obsPoint_502612[0] + 1, obsPoint_502612[1], obsPoint_502612[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_502612[:, 0], obs_measured_502612[:, 1] + botm[nlay-1,int(idx2[1]),int(idx2[2])], color='b', ls=':', label='Pegel 6B (502612)')
plt.plot(obs_measured_502611[:, 0], obs_measured_502611[:, 1] + botm[nlay-1,int(idx2_2[1]),int(idx2_2[2])], color='r', ls=':', label='Pegel 6D (502611)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('time_series_centre.png', dpi=300, bbox_inches='tight')
plt.show()

### Import measured observation point 4
obs_measured_502442  = np.loadtxt('obs_head_time_502442.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_502441  = np.loadtxt('obs_head_time_502441.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_502442 = [plot_layer, 3900+xoff, 5230]
obsPoint_502441 = [0, 3900+xoff, 5230]
## Convert observation point to layer, column, row format
idx4 = (obsPoint_502442[0], 
       round(obsPoint_502442[2]/delr,0), 
       round(obsPoint_502442[1]/delc,0))
idx4_1 = (obsPoint_502441[0], 
       round(obsPoint_502441[2]/delr,0), 
       round(obsPoint_502441[1]/delc,0))

idx_bot = dis.botm.array[:, int(idx4[1]), int(idx4[2])]
ts = headobj.get_ts(idx4)
ts_1 = headobj.get_ts(idx4_1)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = 3900 m and y = 5230 m'.format(obsPoint_502442[0] + 1, obsPoint_502442[1], obsPoint_502442[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_502442[:, 0], obs_measured_502442[:, 1] + botm[nlay-1,int(idx4[1]),int(idx4[2])], ls=':', label='Pegel 6B (502442)')
plt.plot(obs_measured_502441[:, 0], obs_measured_502441[:, 1] + botm[nlay-1,int(idx4_1[1]),int(idx4_1[2])], color="red", ls=':', label='Pegel 6D (502441)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('time_series_south.png', dpi=300, bbox_inches='tight')
plt.show()

### Import measured observation point 5
obs_measured_503762  = np.loadtxt('obs_head_time_503762.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_503761  = np.loadtxt('obs_head_time_503761.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_503762 = [plot_layer, 1400+xoff, 1900]
obsPoint_503761 = [0, 1400+xoff, 1900]
## Convert observation point to layer, column, row format
idx5 = (obsPoint_503762[0], 
       round(obsPoint_503762[2]/delr,0), 
       round(obsPoint_503762[1]/delc,0))
idx5_1 = (obsPoint_503761[0], 
       round(obsPoint_503761[2]/delr,0), 
       round(obsPoint_503761[1]/delc,0))

idx_bot = dis.botm.array[:, int(idx5[1]), int(idx5[2])]
ts = headobj.get_ts(idx5)
ts_1 = headobj.get_ts(idx5_1)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = 3900 m and y = 5230 m'.format(obsPoint_503762[0] + 1, obsPoint_503762[1], obsPoint_503762[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_503762[:, 0], obs_measured_503762[:, 1] + botm[nlay-1,int(idx5[1]),int(idx5[2])], ls=':', label='Pegel 6B (503762)')
plt.plot(obs_measured_503761[:, 0], obs_measured_503761[:, 1] + botm[nlay-1,int(idx4_1[1]),int(idx4_1[2])], color="red", ls=':', label='Pegel 6D (503761)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('time_series_west.png', dpi=300, bbox_inches='tight')
plt.show()

### Import measured observation point 6
obs_measured_805722  = np.loadtxt('obs_head_time_805722.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_805721  = np.loadtxt('obs_head_time_805721.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_805722 = [plot_layer, 1900+xoff, 3400]
obsPoint_805721 = [0, 1900+xoff, 3400]
## Convert observation point to layer, column, row format
idx6 = (obsPoint_805722[0], 
       round(obsPoint_805722[2]/delr,0), 
       round(obsPoint_805722[1]/delc,0))
idx6_1 = (obsPoint_805721[0], 
       round(obsPoint_805721[2]/delr,0), 
       round(obsPoint_805721[1]/delc,0))

idx_bot = dis.botm.array[:, int(idx6[1]), int(idx6[2])]
ts = headobj.get_ts(idx6)
ts_1 = headobj.get_ts(idx6_1)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = 3900 m and y = 5230 m'.format(obsPoint_805722[0] + 1, obsPoint_805722[1], obsPoint_805722[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_805722[:, 0], obs_measured_805722[:, 1] + botm[nlay-1,int(idx6[1]),int(idx6[2])], ls=':', label='Pegel 6B (805722)')
plt.plot(obs_measured_805721[:, 0], obs_measured_805721[:, 1] + botm[nlay-1,int(idx6_1[1]),int(idx6_1[2])], color="red", ls=':', label='Pegel 6D (805721)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
plt.savefig('time_series_southwest.png', dpi=300, bbox_inches='tight')
plt.show()

### Import measured observation point 7
obs_measured_813762  = np.loadtxt('obs_head_time_813762.csv', 
                           delimiter=",",
                           skiprows = 1)
#obs_measured_503761  = np.loadtxt('obs_head_time_503761.csv', 
#                           delimiter=",",
#                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_813762 = [plot_layer, 2100+xoff, 3000]
#obsPoint_503761 = [0, 1400, 1900]
## Convert observation point to layer, column, row format
idx5 = (obsPoint_813762[0], 
       round(obsPoint_813762[2]/delr,0), 
       round(obsPoint_813762[1]/delc,0))
#idx5_1 = (obsPoint_503761[0], 
#       round(obsPoint_503761[2]/delr,0), 
#       round(obsPoint_503761[1]/delc,0))

idx_bot = dis.botm.array[:, int(idx5[1]), int(idx5[2])]
ts = headobj.get_ts(idx5)
#ts_1 = headobj.get_ts(idx5_1)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = 3900 m and y = 5230 m'.format(obsPoint_813762[0] + 1, obsPoint_813762[1], obsPoint_813762[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
#plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
plt.plot(obs_measured_813762[:, 0], obs_measured_813762[:, 1] + botm[nlay-1,int(idx5[1]),int(idx5[2])], ls=':', label='Pegel 6B (813762)')
#plt.plot(obs_measured_503761[:, 0], obs_measured_503761[:, 1] + botm[nlay-1,int(idx4_1[1]),int(idx4_1[2])], color="red", ls=':', label='Pegel 6D (503761)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
#plt.savefig('time_series_west.png', dpi=300, bbox_inches='tight')
plt.show()

### Import Ziellinie 1
obs_target_504943  = np.loadtxt('obs_head_time_504943.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_504943 = [plot_layer, 2529636-xul+xoff, yul-5658103]
## Convert observation point to layer, column, row format
idx = (obsPoint_504943[0], 
       round(obsPoint_504943[2]/delr,0), 
       round(obsPoint_504943[1]/delc,0))

idx_bot = dis.botm.array[:, int(idx[1]), int(idx[2])]
ts = headobj.get_ts(idx)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = 3900 m and y = 5230 m'.format(obsPoint_504943[0] + 1, obsPoint_504943[1], obsPoint_504943[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(obs_target_504943[:, 0], obs_target_504943[:, 1] + botm[nlay-1,int(idx5[1]),int(idx5[2])], ls=':', label='Pegel 6B (504943)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
#plt.savefig('time_series_west.png', dpi=300, bbox_inches='tight')
plt.show()

### Import Ziellinie 2
obs_target_504952  = np.loadtxt('obs_head_time_504952.csv', 
                           delimiter=",",
                           skiprows = 1)
## User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint_504952 = [plot_layer, 2530385-xul+xoff, yul-5658172]
## Convert observation point to layer, column, row format
idx = (obsPoint_504952[0], 
       round(obsPoint_504952[2]/delr,0), 
       round(obsPoint_504952[1]/delc,0))

idx_bot = dis.botm.array[:, int(idx[1]), int(idx[2])]
ts = headobj.get_ts(idx)
plt.subplot(1, 1, 1)
ttl = 'Wasserstand im Modellpunkt x = 3900 m and y = 5230 m'.format(obsPoint_504952[0] + 1, obsPoint_504952[1], obsPoint_504943[2])
plt.title(ttl)
plt.xlabel('Zeit in Tagen')
plt.ylabel('Wasserstand in m')
plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
plt.plot(obs_target_504952[:, 0], obs_target_504952[:, 1] + botm[nlay-1,int(idx5[1]),int(idx5[2])], ls=':', label='Pegel 6B (504952)')
plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.axis([1500, 4500, 0, 180])
#plt.savefig('time_series_west.png', dpi=300, bbox_inches='tight')
plt.show()

##budget graphs
mf_list = flopy.utils.MfListBudget(modelname+".list")
budget = mf_list.get_budget()
#for stress_period in range(2,nper):
#    timestep = nstp[stress_period-1]-1
#    data = mf_list.get_data(kstpkper=(timestep,stress_period))
#    plt.title('water budget after ' + str(stress_period + 1) + ' stress period at ' + str(timestep + 1) + ". timestep\n")
#    plt.bar(data['index'], data['value'])
#    plt.xticks(data['index'], data['name'], rotation=45, size=6)
#    plt.ylabel('m3')
#    plt.show()
#    
budget_incremental, budget_cumulative = mf_list.get_dataframes(start_datetime='31-12-2006')
#plt.plot(budget_incremental["STORAGE_IN"].as_matrix()*perlen, label="In: storage")
#plt.plot(budget_incremental["CONSTANT_HEAD_IN"].as_matrix()*perlen, label="In: constant head")
#plt.plot(budget_incremental["MNW2_IN"].as_matrix()*perlen, label="In: MNW2")
#plt.plot(budget_incremental["STORAGE_OUT"].as_matrix()*perlen, label="Out: storage")
#plt.plot(budget_incremental["CONSTANT_HEAD_OUT"].as_matrix()*perlen, label="Out: constant head")
#plt.plot(budget_incremental["MNW2_OUT"].as_matrix()*perlen, label="Out: MNW2")
bar_width = 0.35
plt.bar(layer3_budget_perStressPeriod['stress_period'], budget_incremental["MNW2_OUT"].as_matrix()*perlen, bar_width, label="QBrunnen")
plt.bar(layer3_budget_perStressPeriod['stress_period'], budget_incremental["CONSTANT_HEAD_OUT"].as_matrix()*perlen, bar_width, color= 'green', label="Out: Rand", bottom=budget_incremental["MNW2_OUT"].as_matrix()*perlen)
plt.bar(layer3_budget_perStressPeriod['stress_period'] + bar_width, budget_incremental["STORAGE_IN"].as_matrix()*perlen, bar_width, color='r', label="QVorrat")
plt.bar(layer3_budget_perStressPeriod['stress_period'] + bar_width, budget_incremental["CONSTANT_HEAD_IN"].as_matrix()*perlen, bar_width, color= 'orange', label="QRand", bottom=budget_incremental["STORAGE_IN"].as_matrix()*perlen)
plt.legend(bbox_to_anchor=(0.4, 0.9), bbox_transform=plt.gcf().transFigure)
plt.title('Wasserbilanz Gesamt (in m3 pro Stressperiode)')
plt.axis([0, 9, 0, 1.8e7])
plt.ylabel('m3')
plt.xlabel('Stressperiode')
plt.xticks(layer3_budget_perStressPeriod['stress_period'])
plt.savefig('budget.png', dpi=300, bbox_inches='tight')
plt.show()

###pump graph
### Import measured pump data
pumping_calculated_6B  = np.loadtxt('pump data.csv', 
                           delimiter=",",
                           skiprows = 1,
                           usecols = (0,3))

#plt.bar(layer3_budget_perStressPeriod['stress_period'], budget_incremental["MNW2_OUT"].as_matrix()*perlen, bar_width, label="Out: Modell")
plt.bar(pumping_calculated_6B [:,0]+bar_width, pumping_calculated_6B [:,1], bar_width, color='r', label="Out: BIOS")
plt.legend(bbox_to_anchor=(0.42, 0.9), bbox_transform=plt.gcf().transFigure)
plt.title('Entnahme aus 6B (in m3 pro Jahr)')
plt.axis([0, 8, 0, 1.8e7])
plt.ylabel('m3')
plt.xlabel('Stressperiode')
plt.xticks(layer3_budget_perStressPeriod['stress_period'])
#plt.savefig('budget.png', dpi=300)
plt.show()

### ModelCrossSection
#fig = plt.figure(figsize=(8, 6))
#ax = fig.add_subplot(1, 1, 1)
#ax.set_title('contour_array() and plot_surface()')
modelxsect = flopy.plot.ModelCrossSection(model=mf, line={'row': 75})
#ct = modelxsect.contour_array(head, masked_values=[999.], head=head, levels=levels)
patches = modelxsect.plot_ibound(head=head)
wt = modelxsect.plot_surface(head, masked_values=[999.], color='blue', lw=1)
linecollection = modelxsect.plot_grid()
plt.title('Profilschnitt in W-O-Richtung mit Grundwasserständen')
plt.ylabel('m')
plt.xlabel('m')
plt.axis([0, 5000, 0, 200])
plt.savefig('xsect.png', dpi=300, bbox_inches='tight')
plt.show()

#t = ax.set_title('Column 6 Cross-Section - Model Grid')

### ModelCrossSection
#fig = plt.figure(figsize=(8, 6))
#ax = fig.add_subplot(1, 1, 1)
#ax.set_title('contour_array() and plot_surface()')
modelysect = flopy.plot.ModelCrossSection(model=mf, line={'column': 70})
#ct = modelxsect.contour_array(head, masked_values=[999.], head=head, levels=levels)
patches = modelysect.plot_ibound(head=head)
wt = modelysect.plot_surface(head, masked_values=[999.], color='blue', lw=1)
linecollection = modelysect.plot_grid()
plt.title('Profilschnitt in N-S-Richtung mit Grundwasserständen')
plt.ylabel('m')
plt.xlabel('m')
plt.axis([0, 5400, 0, 200])
plt.savefig('ysect.png', dpi=300, bbox_inches='tight')
plt.show()
#t = ax.set_title('Column 6 Cross-Section - Model Grid')
    
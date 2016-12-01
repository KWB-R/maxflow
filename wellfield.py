import numpy as np
import flopy
import os
import pandas as pd
import math

# Model domain and grid definition
Lx = 5000.
Ly = 5000.
ztop = 80.
zbot = 0.
nlay = 3
grid_spacing = 50
delr = grid_spacing
delc = grid_spacing
delv = np.array([40, 20, 20], dtype=np.float32)
nrow = int(Lx / delr)
ncol = int(Ly / delc)
botm = np.array([ztop - delv[0],ztop - sum(delv[0:2]), zbot], dtype=np.float32)
hk = np.array([2e-5*3600*24, 10e-8*3600*24, 2e-5*3600*24], #horizontal conductivity
              dtype=np.float32)
vka =  np.array([2e-5*3600*24, 10e-8*3600*24, 2e-5*3600*24], #vertical conductivity
                dtype=np.float32)
sy = np.array([0.123, 0.023, 0.123], #specific yield
              dtype=np.float32)
ss = np.array([1.e-3, 1.e-3, 1.e-3], #specific storage
              dtype=np.float32)
laytyp = np.int_([1, 0, 0]) # 1 - ungespannt, 0 - gespannt

plot_layer = nlay-1

# Variables for the BAS package
# Note that changes from the previous tutorial!
ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
ibound[:, :, 0] = -1
ibound[:, :, -1] = 1


#starting heads
iniHead = 80.
strt = iniHead * np.ones((nlay, nrow, ncol), dtype=np.float32)


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
totsim = 3*365+1


# Time step parameters
nper = 10 #number of stress periods
steady = [True, False, False, False, False, False, False, False, False, False] #type of sress period

t_perPeriod = (totsim-1)/(nper-1) # totsim - 1 (-1 because steady state time period = 1day)
#perlen = [1, 122, 122, 122, 122, 122, 122, 122, 122] #length of a stress period
perlen =np.append([1], np.repeat(t_perPeriod,nper-1)) #length of a stress period
ntimesteps = np.ceil(perlen[1:nper]/max_dt).astype(int)

nstp = np.append([1], ntimesteps) #number of time steps in a stress period

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
                               steady = steady)
                               
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
                               laytyp = laytyp)
pcg = flopy.modflow.ModflowPcg(mf) #Preconditioned Conjugate-Gradient

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

def create_wellfield(layers = [0],
                     spacing_x = 1000, 
                     spacing_y = 1000, 
                     extent_x = [4250,4750],
                     extent_y = [500,4500],
                     pumping_rate = 0):

    layers = np.array(layers)
    extent_x = np.array(extent_x) #m
    extent_y = np.array(extent_y) #m

    rows = np.array(range(extent_y[0],
                      extent_y[1]+1,
                      spacing_y)) / delr
                       
    rows = np.round(rows,0)

    cols = np.array(range(extent_x[0],
                       extent_x[1]+1,
                       spacing_x)) / delc
                       
    cols = np.round(cols,0)
    num_wells = rows.size * cols.size
    wellfield = pd.DataFrame()

    for layer in layers:
        for well_row in rows:
            for well_col in cols:
                tmp = pd.DataFrame([[layer,well_row, well_col, pumping_rate]],
                               columns=['layer','row','column','pumping_rate'])
                wellfield = wellfield.append(tmp,ignore_index=True)
    msg = 'Added ' +  str(num_wells) + ' pumping wells in layer ' + str(layer) + ' each with Q = ' + str(pumping_rate) + ' m3/day' 
    print(msg)
    return(wellfield)
#wellfield_df = pd.DataFrame(wellfield,columns=["layer","row","column","pumping_rate"])              
#pd.DataFrame.as_matrix(wellfield_df)                  


## Well package 
wel_sp1 = create_wellfield()
wel_sp3 = create_wellfield(spacing_x = 1000, spacing_y = 2000, pumping_rate = -2000)
wel_sp4 = create_wellfield(spacing_x =1000, spacing_y = 1000, pumping_rate = -1900)
wel_sp5 = create_wellfield(spacing_x =1000, spacing_y = 1000, pumping_rate = -1600)
wel_sp6 = create_wellfield(spacing_x =1000, spacing_y = 500, pumping_rate = -1200)
wel_sp7 = create_wellfield(spacing_x =500, spacing_y = 250, pumping_rate = -200)
wel_sp8 = create_wellfield(spacing_x =250, spacing_y = 250, pumping_rate = -105)
wel_sp9 = create_wellfield(spacing_x =250, spacing_y = 250, pumping_rate = -90)                                           
wel_sp10 = create_wellfield(spacing_x =250, spacing_y = 250,pumping_rate=-90)                                            
wel_sp11 = create_wellfield(spacing_x =150, spacing_y = 150,pumping_rate=-150)
wel_sp12 = create_wellfield(spacing_x =150, spacing_y = 150,pumping_rate=-70)
                                            
wfs = wel_sp9.append(create_wellfield(spacing_x =1000, spacing_y = 1000, 
                        extent_x = [2500,2600],
                        extent_y = [500,4500],
                        pumping_rate=-50*24))


pumping_data =       {0: pd.DataFrame.as_matrix(wel_sp1), 
                      1: pd.DataFrame.as_matrix(wel_sp1), 
                      2: pd.DataFrame.as_matrix(wel_sp3),
                      3: pd.DataFrame.as_matrix(wel_sp4),
                      4: pd.DataFrame.as_matrix(wel_sp5),
                      5: pd.DataFrame.as_matrix(wel_sp6), 
                      6: pd.DataFrame.as_matrix(wel_sp7),
                      7: pd.DataFrame.as_matrix(wel_sp8),
                      8: pd.DataFrame.as_matrix(wel_sp9),
                      9: pd.DataFrame.as_matrix(wfs)}
#                      
#wel = flopy.modflow.ModflowWel(mf, stress_period_data = pumping_data)

### MNW-2 package

node_data  = pd.read_csv('wells_nodes.csv')

node_data["i"] = (node_data["y"]/delr).astype(int)
node_data["j"] = (node_data["x"]/delc).astype(int)

wells_location = node_data


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
                             
                               
ocdict = {}
for sper in range(0,nper-1):
    if  steady[sper]==False:
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
try:
    mmw2
except NameError:
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
    

# Imports
import matplotlib.pyplot as plt
import flopy.utils.binaryfile as bf

# Get the headfile and budget file objects
headobj = bf.HeadFile(modelname+'.hds')
times = headobj.get_times()
cbb = bf.CellBudgetFile(modelname+'.cbc')

# Setup contour parameters
levels = np.linspace(0, 80, 17)
extent = (delr/2., Lx - delr/2., delc/2., Ly - delc/2.)
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
    pumping_wells = pumping_data[str_per]
    mnw_wells = wells_info[wells_info['per'] == str_per-1] 
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
    try:
        mmw2
    except NameError:
        print("Using MNW2 package and plotting active wells")
        plt.plot(mnw_wells['j']*delc, 
                 mnw_wells['i']*delr, 
                 lw=0, 
                 marker='o', 
                 markersize=3, 
                 markeredgewidth=1,
                 markeredgecolor='black', 
                 markerfacecolor=mfc, 
                 zorder=9)
    else:
        print("Using WEL package and plotting active wells")
        plt.plot(pumping_wells[:,2]*delc, 
                 pumping_wells[:,1]*delr, 
                 lw=0, marker='o', markersize=3, 
                 markeredgewidth=1,
                 markeredgecolor='black', markerfacecolor=mfc, zorder=9)
    plt.show()
    #plt.text(wpt[0]+25, wpt[1]-25, 'well', size=12, zorder=12)


plt.show()



head = headobj.get_data(totim=times[len(times)-1])
levels = np.arange(-50, 10, .5)

il = 0
time = times[len(times)-1]
mytitle = 'Heads in layer ' + str(il) + ' after '+ str(time) + ' days of simulation'
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
title = ax.set_title(mytitle)
modelmap = flopy.plot.ModelMap(model=mf, rotation=0)
quadmesh = modelmap.plot_ibound()
contour_set = modelmap.plot_array(head[il,:,:], 
                                  masked_values=[999.], 
                                  alpha=0.5)
linecollection = modelmap.plot_grid()
cb = plt.colorbar(contour_set, shrink=0.4)
plt.show()


###Plot the head versus time

### Import measured observation point

obs_measured  = np.loadtxt('obs_head_time.csv', 
                           delimiter=",",
                           skiprows = 1)

### User defined observation point in format: layer, x, y-coordinate (absolute)
obsPoint = [plot_layer, 4850, 2500]

#### Convert observation point to layer, column, row format
idx = (obsPoint[0], 
       round(obsPoint[2]/delr,0), 
       round(obsPoint[1]/delc,0))
ts = headobj.get_ts(idx)
plt.subplot(1, 1, 1)
ttl = 'Head in layer {0} at x = {1} m and y = {2} m'.format(obsPoint[0] + 1, obsPoint[1], obsPoint[2])
plt.title(ttl)
plt.xlabel('time in days')
plt.ylabel('head in m')
plt.plot(ts[:, 0], ts[:, 1], label='simulated')
plt.plot(obs_measured[:, 0], obs_measured[:, 1], color="red", label='measured (814193)')
plt.legend()
plt.show()

mf_list = flopy.utils.MfListBudget(modelname+".list")
budget = mf_list.get_budget()
for stress_period in range(2,nper):
    timestep = nstp[stress_period-1]-1
    data = mf_list.get_data(kstpkper=(timestep,stress_period))
    plt.title('water budget for ' + str(stress_period + 1) + ' stress period at ' + str(timestep + 1) + ". timestep\n")
    plt.bar(data['index'], data['value'])
    plt.xticks(data['index'], data['name'], rotation=45, size=6)
    plt.ylabel('m3')
    plt.show()
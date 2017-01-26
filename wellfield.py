import numpy as np
import os
import flopy
import pandas as pd
import flopy.utils.binaryfile as bf
from create_model import *
from analyse_model import *
from get_layerBudget import *

# 0) Creating MNW file with R (still open to be integrated!)

Ly = 5400.
Lx = 4000. 

upleft_coord = calc_model_wellcoordinates(Ly = Ly, 
                                          Lx = Lx, 
                                          csvDir = '.',
                                          csvFile = 'wells_nodes.csv',
                                          exp_dir = '.')

xul = upleft_coord['xul']
yul = upleft_coord['yul']
proj4_str = 'ESPG:31466'

###############################################################################
##### Scenario 1: Combined leakage with constant kf for 9 stress periods
###############################################################################


print("###############################################################################")
print("Scenario 1: Combined leakage with constant kf for 9 stress periods")
print("###############################################################################\n")
# 1) Create model 
print("Step 1: Create MODFLOW model")
mf = create_model(Ly = Ly, 
                  Lx = Lx,
                  ztop = 200.,
             zbot_north = 0.,
             nlay = 3,
             grid_spacing = 50,
             delv = np.array([160, 20, 20], dtype = np.float32),
             botm_gradient_northSouth = 65/5400, 
             head_north = [160, 160, 110],
             head_gradient_northSouth = -40/5400, 
             hk = np.array([2e-5*3600*24, 1e-9*3600*24, 3e-5*3600*24], #horizontal conductivity
                           dtype=np.float32),
             vka =  np.array([2e-5*3600*24, 1e-9*3600*24, 3e-5*3600*24], #vertical conductivity
                                 dtype=np.float32),
             area_borehole = 0.3, ###meter^2
             kf_borehole = 1e-4*24*3600, #### meter / day
             sy = np.array([0.123, 0.023, 0.123], #specific yield
                               dtype=np.float32),
             ss = np.array([1.e-4, 1.e-4, 1.e-4], #specific storage
                               dtype=np.float32),
             laytyp = np.int_([1, 1, 1]), # 1 - ungespannt, 0 - gespannt):
             totsim = 9*365, ### Desired total simulation time
             nper = 9, #number of stress periods
             hclose = 1E-4, 
             rclose = 5E-4, 
             constantcv = True,
             each_time_step = True, ### if True (for all time steps), if False only for last time 
                 #for each stress period 
             modelname = 'wellfield', 
             model_ws = '.',
             xul = xul,
             yul = yul,
             proj4_str = proj4_str,
             start_datetime = '1/1/2007')
                 
                 
# 2) Run the model
print("Step 2: Run MODFLOW model")
success, mfoutput = mf.run_model(silent=False, pause=False)
if not success:
    raise Exception('MODFLOW did not terminate normally.')
    
# 3) Analyse model results (see analyse_model.py)
print("Step 3: Analyse MODFLOW results (see analyse_model.py)")
analyse_model(modelname = 'wellfield', 
              plot_layer = 2,
              obs_start_date = '1/01/2007')

                 

###############################################################################
##### Scenario 2: Combined leakage with varying kf for each of the 9 stress periods
#####             (dep)
###############################################################################   
print("###############################################################################")
print("Scenario 2: Combined leakage with varying kf for each of the 9 stress periods")
print("###############################################################################\n")
# 0) Creating well_times.csv & well_nodes.csv for each stress period 
### (based on well_times.csv & well_nodes.csv in project folder root)

create_mnw2_csv_perPeriod(csvdir = '.',  ### location of org.  well_times.csv & well_nodes.csv (default: project root)
                          basedir = 'SP' ### base subdir for exporting both files for each stress period (default: 'SP')
                          )

#### Define constant model parameters:
ztop = 200.
zbot_north = 0.
nlay = 3
grid_spacing = 50
delv = np.array([160, 20, 20], dtype = np.float32)
botm_gradient_northSouth = 65/5400
hk = np.array([2e-5*3600*24, 1e-9*3600*24, 3e-5*3600*24], #horizontal conductivity
                       dtype=np.float32)
vka =  np.array([2e-5*3600*24, 1e-9*3600*24, 3e-5*3600*24], #vertical conductivity
                             dtype=np.float32)
area_borehole = 0.3 ###meter^2
kf_borehole = 1e-4*24*3600 #### meter / day
sy = np.array([0.123, 0.023, 0.123], #specific yield
               dtype=np.float32)
ss = np.array([1.e-4, 1.e-4, 1.e-4], #specific storage
               dtype=np.float32)
laytyp = np.int_([1, 1, 1]) # 1 - ungespannt, 0 - gespannt):
totsim = 1*365 ### Desired total simulation time
nper = 1 #number of stress periods
hclose = 1E-4
rclose = 5E-4
constantcv = True
each_time_step = True ### if True (for all time steps), if False only for last time for each stress period 
modelname = 'wellfield'


for stress_period in range(0,9):
    ### Define model folder for each stress period
    model_ws = 'SP' + str(stress_period)
    start_datetime = '1/1/' + str(2007 + stress_period)
    # 1) Create model 
    print("Step 1: Create MODFLOW model for stress period " + str(stress_period))

    if stress_period == 0:
        mf = create_model(Ly = Ly,
             Lx = Lx, 
             ztop = ztop,
             zbot_north = zbot_north,
             nlay = nlay,
             grid_spacing = grid_spacing,
             delv = delv,
             botm_gradient_northSouth = 65/5400, 
             head_north = [160, 160, 110], ### will be used if 'head_array' is 'None'  
             head_gradient_northSouth = -40/5400, 
             head_array = None, ### takes an array as starting hea, if 'None: head_north and/or head_gradient_northSouth will be used
             hk = hk,
             vka =  vka,
             area_borehole = area_borehole, 
             kf_borehole =  kf_borehole, 
             sy = sy,
             ss = ss, 
             laytyp = laytyp, 
             totsim = totsim, 
             nper = nper,
             hclose = hclose, 
             rclose = rclose, 
             constantcv = constantcv,
             each_time_step =  each_time_step, 
             modelname = modelname, 
             model_ws = model_ws,
             xul = xul,
             yul = yul,
             proj4_str = proj4_str, 
             start_datetime = start_datetime)
                 
    
    else:
        old_model_dir = os.path.dirname(mf.dis.fn_path)
        headobj = bf.HeadFile(os.path.join(old_model_dir, modelname+'.hds'))
        times = headobj.get_times()
        head = headobj.get_data(totim=times[len(times)-1])
        mf = create_model(Ly = Ly,
             Lx = Lx, 
             ztop = ztop,
             zbot_north = zbot_north,
             nlay = nlay,
             grid_spacing = grid_spacing,
             delv = delv,
             botm_gradient_northSouth = 65/5400, 
             head_north = None, ### will be used if 'head_array' is 'None'  
             head_gradient_northSouth = -40/5400, 
             head_array = head, ### takes an array as starting hea, if 'None: head_north and/or head_gradient_northSouth will be used
             hk = hk,
             vka =  vka,
             area_borehole = area_borehole, 
             kf_borehole =  kf_borehole, 
             sy = sy,
             ss = ss, 
             laytyp = laytyp, 
             totsim = totsim, 
             nper = nper,
             hclose = hclose, 
             rclose = rclose, 
             constantcv = constantcv,
             each_time_step =  each_time_step, 
             modelname = modelname, 
             model_ws = model_ws,
             xul = xul,
             yul = yul,
             proj4_str = proj4_str, 
             start_datetime = start_datetime)
                 
    # 2) Run the model
    print("Step 2: Run MODFLOW model for stress period " + str(stress_period))
    success, mfoutput = mf.run_model(silent=False, pause=False)
    if not success:
        raise Exception('MODFLOW did not terminate normally.')
    
    # 3) Analyse model results (see analyse_model.py)
    print("Step 3: Analyse MODFLOW results (see analyse_model.py)")
    analyse_model(modelname = modelname, 
                  model_ws = model_ws,
                  plot_layer = 2,
                  obs_start_date = '1/01/2007')


            



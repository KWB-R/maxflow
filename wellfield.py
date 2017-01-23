import numpy as np
import flopy
import pandas as pd
from create_model import *
from analyse_model import *
from get_layerBudget import *

# 0) Creating MNW file with R (still open to be integrated!)
#import subprocess
#
## Define command and arguments
#command = 'C:\\Program Files\\R\\R-3.3.2\\bin\\Rscript.exe'
#path2script = 'C:\\Users\\mrustl\\Desktop\\WC_Maxflow\\branches\\code_refactoring\\create_mnw2_csvs.R'
#
## Variable number of args in a list
#args = ['5400', '4000', 'FALSE', 'FALSE']
#
## Build subprocess command
#cmd = [command, path2script] + args
#
## check_output will run the command and store to result
#x = subprocess.check_output(cmd, universal_newlines=True)



# 1) Create model 
print("Step 1: Create MODFLOW model")
mf = create_model(Ly = 5400.,
             Lx = 4000., 
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
             modelname = 'wellfield')



# 2) Run the model
print("Step 2: Run MODFLOW model")
success, mfoutput = mf.run_model(silent=False, pause=False)
if not success:
    raise Exception('MODFLOW did not terminate normally.')

# 3) Analyse model results (see analyse_model.py)
print("Step 3: Analyse MODFLOW results (see analyse_model.py)")
analyse_model(modelname = 'wellfield', 
              plot_layer = 2)

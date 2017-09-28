import os
import numpy as np
import flopy
import pandas as pd
import sys
from create_model import calc_model_wellcoordinates

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

### 81419
xkoord_1= 2529918
ykoord_1= 5661636
obs_measured_1_6B  = np.loadtxt('obs_head_time_814193.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_1_6D  = np.loadtxt('obs_head_time_814192.csv', 
                           delimiter=",",
                           skiprows = 1)

###80572
xkoord_2= 2530304
ykoord_2= 5658699
obs_measured_2_6B  = np.loadtxt('obs_head_time_805722.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_2_6D  = np.loadtxt('obs_head_time_805721.csv', 
                           delimiter=",",
                           skiprows = 1)


###50495
xkoord_3=  2530385
ykoord_3=  5658172
obs_target_3_6B  = np.loadtxt('target_head_time_504952.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_3_6B  = np.loadtxt('obs_head_time_504952.csv', 
                           delimiter=",",
                           skiprows = 1)
#obs_measured_3_6D  = np.loadtxt('obs_head_time_504953.csv', 
#                           delimiter=",",
#                           skiprows = 1)
###50376
xkoord_4=  2529833
ykoord_4=  5660219
obs_measured_4_6B  = np.loadtxt('obs_head_time_503762.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_4_6D  = np.loadtxt('obs_head_time_503761.csv', 
                           delimiter=",",
                           skiprows = 1)

###50337
xkoord_5= 2530743
ykoord_5= 5661714
obs_measured_5_6B  = np.loadtxt('obs_head_time_503372.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_5_6D  = np.loadtxt('obs_head_time_503371.csv', 
                           delimiter=",",
                           skiprows = 1)

###50261
xkoord_6=  2530515
ykoord_6=  5659857
obs_measured_6_6B  = np.loadtxt('obs_head_time_502612.csv', 
                           delimiter=",",
                           skiprows = 1)
obs_measured_6_6D  = np.loadtxt('obs_head_time_502611.csv', 
                           delimiter=",",
                           skiprows = 1)


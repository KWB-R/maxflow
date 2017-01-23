# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 10:58:55 2017

@author: mrustl
"""

import numpy as np
import pandas as pd
import flopy.utils.binaryfile as bf



def get_layerbudget(modelname,
                    nper, 
                    perlen,
                    nlay, 
                    debug = True
                    ):
    bud_agg = pd.DataFrame(columns=['stress_period',
                              'time_step',
                              'layer',  
                              'STORAGE_IN',
                              'STORAGE_OUT',
                              'CONSTANT_HEAD_IN',
                              'CONSTANT_HEAD_OUT',
                              'FLOW_RIGHT_FACE', 
                              'FLOW_FRONT_FACE',
                              'FLOW_LOWER_FACE'])
    cbb = bf.CellBudgetFile(modelname+'.cbc')    
    for stress_period in np.arange(0,nper).astype('int'):
        for time_step in np.arange(0,perlen[stress_period]).astype('int'):
            bud = cbb.get_data(kstpkper = (time_step,stress_period), 
                       full3D=True)    
            for layer in np.arange(0,nlay).astype('int'): 
                if debug: print("Stress period: " + str(stress_period + 1) + ", Time step: " + str(time_step + 1) + ", Layer: " + str(layer + 1))
                tmp = pd.DataFrame([[stress_period, 
                         time_step + 1, 
                         layer,
                         bud[0][layer][bud[0][layer]>0].sum(),
                         bud[0][layer][bud[0][layer]<0].sum(),
                         bud[1][layer][bud[1][layer]>0].sum(),
                         bud[1][layer][bud[1][layer]<0].sum(),
                         bud[2][layer].sum(),
                         bud[3][layer].sum(),
                         bud[4][layer].sum()
                         ]],
                     columns=['stress_period',
                              'time_step',
                              'layer',  
                              'STORAGE_IN',
                              'STORAGE_OUT',
                              'CONSTANT_HEAD_IN',
                              'CONSTANT_HEAD_OUT',
                              'FLOW_RIGHT_FACE', 
                              'FLOW_FRONT_FACE',
                              'FLOW_LOWER_FACE'])
                bud_agg = bud_agg.append(tmp,ignore_index=True)
    bud_agg.loc[:,['CONSTANT_HEAD_IN']] = bud_agg['CONSTANT_HEAD_IN'].as_matrix().astype("float32")
    bud_agg.loc[:,['CONSTANT_HEAD_OUT']] = bud_agg['CONSTANT_HEAD_OUT'].as_matrix().astype("float32")
    bud_agg.loc[np.isnan(bud_agg['CONSTANT_HEAD_IN']),['CONSTANT_HEAD_IN']] = 0
    bud_agg.loc[np.isnan(bud_agg['CONSTANT_HEAD_OUT']),['CONSTANT_HEAD_OUT']] = 0
    return(bud_agg)


# Imports
import matplotlib.pyplot as plt
import os
import flopy
import flopy.utils.binaryfile as bf
from get_layerBudget import *


def getStressPeriod(time, 
                    stress_period_ende):
        
    cum_perende = np.cumsum(stress_period_ende)
        
    stress_period = sum(time > cum_perende)
       
    return(stress_period )
    
    
def get_mnw2_info(mf):
    delc =  int(min(mf.dis.delc.array))
    delr = int(min(mf.dis.delr.array))
    node_data  = pd.read_csv('wells_nodes.csv')

    node_data["k"] = node_data["k"].astype(int)
    node_data["j"] = (node_data["x"]/delc).astype(int)
    node_data["i"] = (node_data["y"]/delr).astype(int)
    node_data["j"] = (node_data["x"]/delc).astype(int)

    wells_location = node_data

    stress_period_data  = pd.read_csv('wells_times.csv')

    wells_info =  pd.merge(left=wells_location,
                       right=stress_period_data, 
                       left_on='wellid', 
                        right_on='wellid')
    wells_info = wells_info[wells_info['qdes'] != 0]

    return(wells_info)
        

def analyse_model(modelname = 'wellfield',
                  plot_layer = 2):

    m = flopy.modflow.Modflow.load(modelname + '.nam')
    
       
    wells_info = get_mnw2_info(m)
    
     # Get the headfile and budget file objects
    headobj = bf.HeadFile(modelname+'.hds')
    times = headobj.get_times()
    cbb = bf.CellBudgetFile(modelname+'.cbc')

    nlay = m.dis.nlay
    nper = m.dis.nper
    perlen = m.dis.perlen.array
    nstp = m.dis.nstp.array
    ### Get layer based budget for each time step
    layer_budget = get_layerbudget(modelname, 
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
    plt.bar(layer3_budget_perStressPeriod['stress_period'], layer3_budget_perStressPeriod['MNW2_OUT'], bar_width, color='b', label="QBrunnen")
    plt.bar(layer3_budget_perStressPeriod['stress_period'], 
        layer3_budget_perStressPeriod['CONSTANT_HEAD_OUT'], 
        bar_width, 
        color= 'green', 
        label="Out: Rand", 
        bottom=layer3_budget_perStressPeriod['MNW2_OUT'])
    plt.bar(layer3_budget_perStressPeriod['stress_period'], 
        layer3_budget_perStressPeriod['STORAGE_OUT'], 
        bar_width, 
        color= 'brown', 
        label="Out: Vorrat", 
        bottom=layer3_budget_perStressPeriod['MNW2_OUT']+layer3_budget_perStressPeriod['CONSTANT_HEAD_OUT'])
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
    plt.legend(bbox_to_anchor=(0.4, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.title('Entnahme aus 6B (in m3 pro Stressperiode)')
    plt.axis([0, 10, 0, 1.8e7])
    plt.ylabel('m3')
    plt.xlabel('Stressperiode')
    plt.xticks(layer3_budget_perStressPeriod['stress_period'])
    plt.savefig('budget_layer3.png', dpi=300, bbox_inches='tight')
    plt.show()


    # Setup contour parameters
    levels = np.linspace(0, 80, 17)  ### dummy values need to be modified!!!!
    
    delr = int(m.dis.delr.array[0])  
    delc = int(m.dis.delc.array[0])
    Ly = delr * m.dis.nrow
    Lx = delc * m.dis.ncol
    
    extent = (delr/2., Ly - delr/2., delc/2., Lx - delc/2.)
    print('Levels: ', levels)
    print('Extent: ', extent)
    #
    ## Well point
    #
    #def getStressPeriod(time, 
    #                    stress_period_ende = perlen):
    #    
    #    cum_perende = np.cumsum(stress_period_ende)
    #    
    #    stress_period = sum(time > cum_perende)
    #   
    #    return(stress_period )
    #    
    #
    ## Make the plots
    #for iplot, time in enumerate(times):
    #    print('*****Processing time: ', time)
    #    
    #    ### Get head data
    #    head = headobj.get_data(totim=time)
    #    
    #    ### Get pumping wells 
    #    str_per = getStressPeriod(time)
    # 
    #    #Print statistics
    #    print('Head statistics')
    #    print('  min: ', head.min())
    #    print('  max: ', head.max())
    #    print('  std: ', head.std())
    #
    #    # Extract flow right face and flow front face
    #    frf = cbb.get_data(text='FLOW RIGHT FACE', totim=time)[0]
    #    fff = cbb.get_data(text='FLOW FRONT FACE', totim=time)[0]
    #
    #    #Create the plot
    #    #plt.subplot(1, len(mytimes), iplot + 1, aspect='equal')
    #    plt.subplot(1, 1, 1, aspect='equal')
    #    plt.title('total time: ' + str(time) + ' days\n(Stress period: ' + str(str_per) + ")")
    #    modelmap = flopy.plot.ModelMap(model=m, layer=plot_layer)
    #    qm = modelmap.plot_ibound()
    #    lc = modelmap.plot_grid()
    #    #qm = modelmap.plot_bc('GHB', alpha=0.5)
    #    cs = modelmap.contour_array(head, levels=levels)
    #    plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f', zorder=11)
    #    quiver = modelmap.plot_discharge(frf, fff, head=head)
    ##    mfc = 'None'
    ##    if (iplot+1) == len(times):
    #    mfc='black'
    #    plt.plot(2050,4950, 
    #                 lw=0, 
    #                 marker='o', 
    #                 markersize=3, 
    #                 markeredgewidth=1,
    #                 markeredgecolor='red', 
    #                 markerfacecolor=mfc, 
    #                 zorder=9)
    #    if 'mnw2' in locals():
    #        print("Using MNW2 package and plotting active wells")
    #        mnw_wells = wells_info[wells_info['per'] == str_per] 
    #        plt.plot(mnw_wells['j']*delc, 
    #                 Ly-(mnw_wells['i']*delr), 
    #                 lw=0, 
    #                 marker='o', 
    #                 markersize=3, 
    #                 markeredgewidth=1,
    #                 markeredgecolor='black', 
    #                 markerfacecolor=mfc, 
    #                 zorder=9)
    #    else:
    #        print("Using WEL package and plotting active wells")
    #        pumping_wells = pumping_data[str_per]
    #        plt.plot(pumping_wells[:,2]*delc, 
    #                 pumping_wells[:,1]*delr, 
    #                 lw=0, marker='o', markersize=3, 
    #                 markeredgewidth=1,
    #                 markeredgecolor='black', markerfacecolor=mfc, zorder=9)
    #    plt.show()
    #    #plt.text(wpt[0]+25, wpt[1]-25, 'well', size=12, zorder=12)
    #
    #
    #plt.show()
    
    
    
    head = headobj.get_data(totim=times[len(times)-1])
    
    time = times[len(times)-1]
    mytitle = 'Heads in layer ' + str(plot_layer) + ' after '+ str(time) + ' days of simulation'
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    title = ax.set_title(mytitle)
    modelmap = flopy.plot.ModelMap(model=m, rotation=0)
    quadmesh = modelmap.plot_ibound()
    contour_set = modelmap.plot_array(head[plot_layer,:,:], 
                                      masked_values=[-1e+30], 
                                      alpha=0.5)
    cs = modelmap.contour_array(head, levels=levels)
    plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f', zorder=11)
    linecollection = modelmap.plot_grid()
    cb = plt.colorbar(contour_set, shrink=0.4)
    mfc = 'None'
    plt.plot(2050,4950, 
                     lw=0, 
                     marker='o', 
                     markersize=3, 
                     markeredgewidth=1,
                     markeredgecolor='red', 
                     markerfacecolor=mfc, 
                     zorder=9)
    plt.plot(2650,3250, 
                     lw=0, 
                     marker='o', 
                     markersize=3, 
                     markeredgewidth=1,
                     markeredgecolor='red', 
                     markerfacecolor=mfc, 
                     zorder=9)
    plt.plot(3900,270, 
                     lw=0, 
                     marker='o', 
                     markersize=3, 
                     markeredgewidth=1,
                     markeredgecolor='red', 
                     markerfacecolor=mfc, 
                     zorder=9)
    plt.plot(2875,5000, 
                     lw=0, 
                     marker='o', 
                     markersize=3, 
                     markeredgewidth=1,
                     markeredgecolor='red', 
                     markerfacecolor=mfc, 
                     zorder=9)
    print("Using MNW2 package and plotting active wells")
    sper = getStressPeriod(time, perlen)
    mnw_wells = wells_info[wells_info['per'] == sper] 
    plt.plot(mnw_wells['j']*delc, 
             Ly-(mnw_wells['i']*delr), 
             lw=0, 
             marker='o', 
             markersize=3, 
             markeredgewidth=1,
             markeredgecolor='black', 
             markerfacecolor=mfc, 
             zorder=9)
    plt.show()
    
    ##Restwasser
    botm =  m.dis.botm.array
    mytitle = 'Restwassermächtigkeit im 6B in 2015' # + str(plot_layer) + ' bottom elevation after '+ str(time) + ' days of simulation'
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, aspect='equal')
    title = ax.set_title(mytitle)
    modelmap = flopy.plot.ModelMap(model=m, rotation=0)
    quadmesh = modelmap.plot_ibound()
    contour_set = modelmap.plot_array(head[plot_layer,:,:]-botm[plot_layer], 
                                      masked_values=[-1e+30], 
                                      alpha=0.5)
    #levels = np.linspace(0, 80, 17)
    #extent = (delr/2., Ly - delr/2., delc/2., Lx - delc/2.)
    #print('Levels: ', levels)
    #print('Extent: ', extent)
    #cs = modelmap.contour_array(head[plot_layer,:,:]-botm[plot_layer], levels=levels)
    #plt.clabel(cs, inline=1, fontsize=10, fmt='%1.1f', zorder=11)
    linecollection = modelmap.plot_grid()
    cb = plt.colorbar(contour_set, shrink=0.4)
    plt.plot(2050,4950, 
                     lw=0, 
                     marker='o', 
                     markersize=8, 
                     markeredgewidth=1,
                     markeredgecolor='red', 
                     markerfacecolor='red', 
                     zorder=9)
    plt.plot(2650,3250, 
                     lw=0, 
                     marker='o', 
                     markersize=8, 
                     markeredgewidth=1,
                     markeredgecolor='red', 
                     markerfacecolor='red', 
                     zorder=9)
    plt.plot(3900,270, 
                     lw=0, 
                     marker='o', 
                     markersize=8, 
                     markeredgewidth=1,
                     markeredgecolor='red', 
                     markerfacecolor='red', 
                     zorder=9)
    plt.plot(2875,5000, 
                     lw=0, 
                     marker='o', 
                     markersize=8, 
                     markeredgewidth=1,
                     markeredgecolor='red', 
                     markerfacecolor='red', 
                     zorder=9)
#    if 'mnw2' in locals():
    print("Using MNW2 package and plotting active wells")
    sper = getStressPeriod(time, m.dis.perlen.array)
    mnw_wells = wells_info[wells_info['per'] == sper] 
    plt.plot(mnw_wells['j']*delc, 
             Ly-(mnw_wells['i']*delr), 
             lw=0, 
             marker='o', 
             markersize=3, 
             markeredgewidth=1,
             markeredgecolor='black', 
             markerfacecolor=mfc, 
             zorder=9)
    plt.savefig('Restwasser.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    ###Plot the head versus time
    
    ### Import measured observation point 1
    obs_measured1  = np.loadtxt('obs_head_time_814193.csv', 
                               delimiter=",",
                               skiprows = 1)
    obs_measured1_1  = np.loadtxt('obs_head_time_814192.csv', 
                               delimiter=",",
                               skiprows = 1)
    ## User defined observation point in format: layer, x, y-coordinate (absolute)
    obsPoint1 = [plot_layer, 2050, 450]
    obsPoint1_1 = [0, 2050, 450]
    
    ## Convert observation point to layer, column, row format
    idx1 = (obsPoint1[0], 
           round(obsPoint1[2]/delr,0), 
           round(obsPoint1[1]/delc,0))
    idx1_1 = (obsPoint1_1[0], 
           round(obsPoint1_1[2]/delr,0), 
           round(obsPoint1_1[1]/delc,0))
    
    idx_bot = m.dis.botm.array[:, int(idx1[1]), int(idx1[2])]
    
    ts = headobj.get_ts(idx1) 
    ts_1 = headobj.get_ts(idx1_1) 
    plt.subplot(1, 1, 1)
    ttl = 'Wasserstand im Modellpunkt x = {1} m and y = {2} m'.format(obsPoint1[0] + 1, obsPoint1[1], obsPoint1[2])
    plt.title(ttl)
    plt.xlabel('Zeit in Tagen')
    plt.ylabel('Wasserstand in m')
    plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
    plt.plot(ts_1[:, 0], ts_1[:, 1], color="red", label='Modell 6D')
    plt.plot(obs_measured1[:, 0], obs_measured1[:, 1] + botm[nlay-1,int(idx1[1]),int(idx1[2])], color="blue", ls=':', label='Pegel 6B (814193)')
    plt.plot(obs_measured1_1[:, 0], obs_measured1_1[:, 1] + botm[nlay-1,int(idx1_1[1]),int(idx1_1[2])], color="red", ls=':', label='Pegel 6D (814192)')
    plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
    plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
    plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.axis([0, 3500, 0, 160])
    plt.savefig('time_series_north.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    ### Import measured observation point 2
    obs_measured_503372  = np.loadtxt('obs_head_time_503372.csv', 
                               delimiter=",",
                               skiprows = 1)
    
    ## User defined observation point in format: layer, x, y-coordinate (absolute)
    obsPoint_503372 = [plot_layer, 2875, 400]
    
    ## Convert observation point to layer, column, row format
    idx2 = (obsPoint_503372[0], 
           round(obsPoint_503372[2]/delr,0), 
           round(obsPoint_503372[1]/delc,0))
    
    idx_bot = m.dis.botm.array[:, int(idx2[1]), int(idx2[2])]
                             
    ts = headobj.get_ts(idx2) 
    plt.subplot(1, 1, 1)
    ttl = 'Wasserstand im Modellpunkt x = 2875 m and y = 400 m'.format(obsPoint1[0] + 1, obsPoint1[1], obsPoint1[2])
    plt.title(ttl)
    plt.xlabel('Zeit in Tagen')
    plt.ylabel('Wasserstand in m')
    plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
    plt.plot(obs_measured_503372[:, 0], obs_measured_503372[:, 1] + botm[nlay-1,int(idx2[1]),int(idx2[2])], ls=':', label='Pegel 6B (503372)')
    plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
    plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
    plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
    plt.legend()
    plt.axis([0, 3500, 0, 160])
    plt.savefig('time_series_east.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    ### Import measured observation point 3
    obs_measured_502612  = np.loadtxt('obs_head_time_502612.csv', 
                               delimiter=",",
                               skiprows = 1)
    
    ## User defined observation point in format: layer, x, y-coordinate (absolute)
    obsPoint_502612 = [plot_layer, 2650, 2230]
    
    ## Convert observation point to layer, column, row format
    idx2 = (obsPoint_502612[0], 
           round(obsPoint_502612[2]/delr,0), 
           round(obsPoint_502612[1]/delc,0))
    
    idx_bot = m.dis.botm.array[:, int(idx2[1]), int(idx2[2])]
    ts = headobj.get_ts(idx2) 
    plt.subplot(1, 1, 1)
    ttl = 'Wasserstand im Modellpunkt x = 2650 m and y = 2230 m'.format(obsPoint1[0] + 1, obsPoint1[1], obsPoint1[2])
    plt.title(ttl)
    plt.xlabel('Zeit in Tagen')
    plt.ylabel('Wasserstand in m')
    plt.plot(ts[:, 0], ts[:, 1], color="blue", label='Modell 6B')
    plt.plot(obs_measured_502612[:, 0], obs_measured_502612[:, 1] + botm[nlay-1,int(idx2[1]),int(idx2[2])], ls=':', label='Pegel 6B (502612)')
    plt.axhline(y=idx_bot[0], color='grey', linestyle='-')
    plt.axhline(y=idx_bot[1], color='grey', linestyle='-')
    plt.axhline(y=idx_bot[2], color='grey', linestyle='-')
    plt.legend()
    plt.axis([0, 3500, 0, 160])
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
    obsPoint_502442 = [plot_layer, 3900, 5230]
    obsPoint_502441 = [0, 3900, 5230]
    ## Convert observation point to layer, column, row format
    idx4 = (obsPoint_502442[0], 
           round(obsPoint_502442[2]/delr,0), 
           round(obsPoint_502442[1]/delc,0))
    idx4_1 = (obsPoint_502441[0], 
           round(obsPoint_502441[2]/delr,0), 
           round(obsPoint_502441[1]/delc,0))
    
    idx_bot = m.dis.botm.array[:, int(idx4[1]), int(idx4[2])]
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
    plt.axis([0, 3500, 0, 160])
    plt.savefig('time_series_south.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    #budget graphs
    mf_list = flopy.utils.MfListBudget(modelname+".list")
    budget = mf_list.get_budget()
    for stress_period in range(0,nper):
        print(stress_period)
        timestep = nstp[stress_period]-1
        data = mf_list.get_data(kstpkper=(timestep,stress_period))
        plt.title('water budget after ' + str(stress_period + 1) + ' stress period at ' + str(timestep + 1) + ". timestep\n")
        plt.bar(data['index'], data['value'])
        plt.xticks(data['index'], data['name'], rotation=45, size=6)
        plt.ylabel('m3')
        plt.show()
    
    perlen = m.dis.perlen.array
    stpers = np.arange(0,nper)
    
    ###pump graph
    ### Import measured pump data
    pumping_calculated_6B  = np.loadtxt('pump data.csv', 
                               delimiter=",",
                               skiprows = 1,
                               usecols = (0,3))
    
    pumping_modelled = np.append(budget_cumulative['MNW2_OUT'][0], 
                                  budget_cumulative['MNW2_OUT'].diff().as_matrix()[1:])
    
    plt.bar(stpers,  pumping_modelled, bar_width, label="Out: Modell")
    plt.bar(pumping_calculated_6B [:,0]+bar_width, pumping_calculated_6B [:,1], bar_width, color='r', label="Out: BIOS")
    plt.legend(bbox_to_anchor=(0.42, 0.9), bbox_transform=plt.gcf().transFigure)
    plt.title('Entnahme aus 6B (in m3 pro Jahr)')
    plt.axis([0, nper, 0, 1.8e7])
    plt.ylabel('m3')
    plt.xlabel('Stressperiode')
    plt.xticks(stpers)
    #plt.savefig('budget.png', dpi=300)
    plt.show()
    
    ztop = max(m.dis.top.array[0])
    ### ModelCrossSection
    #fig = plt.figure(figsize=(8, 6))
    #ax = fig.add_subplot(1, 1, 1)
    #ax.set_title('contour_array() and plot_surface()')
    modelxsect = flopy.plot.ModelCrossSection(model=m, line={'row': 5})
    #ct = modelxsect.contour_array(head, masked_values=[999.], head=head, levels=levels)
    patches = modelxsect.plot_ibound(head=head)
    wt = modelxsect.plot_surface(head, masked_values=[999.], color='blue', lw=1)
    linecollection = modelxsect.plot_grid()
    plt.title('Profilschnitt in W-O-Richtung mit Grundwasserständen')
    plt.ylabel('m')
    plt.xlabel('m')
    plt.axis([0, Lx, 0, ztop])
    plt.savefig('xsect.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    #t = ax.set_title('Column 6 Cross-Section - Model Grid')
    
    ### ModelCrossSection
    #fig = plt.figure(figsize=(8, 6))
    #ax = fig.add_subplot(1, 1, 1)
    #ax.set_title('contour_array() and plot_surface()')
    modelysect = flopy.plot.ModelCrossSection(model=m, line={'column': 70})
    #ct = modelxsect.contour_array(head, masked_values=[999.], head=head, levels=levels)
    patches = modelysect.plot_ibound(head=head)
    wt = modelysect.plot_surface(head, masked_values=[999.], color='blue', lw=1)
    linecollection = modelysect.plot_grid()
    plt.title('Profilschnitt in N-S-Richtung mit Grundwasserständen')
    plt.ylabel('m')
    plt.xlabel('m')
    plt.axis([0, Ly, 0, ztop])
    plt.savefig('ysect.png', dpi=300, bbox_inches='tight')
    plt.show()
    #t = ax.set_title('Column 6 Cross-Section - Model Grid')
    
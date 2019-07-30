# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 16:09:26 2018

Plotting scripts for Plant_Designer

@author: Drizdar
"""

import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import polyfit
import numpy as np
import pro_sys_eqns as pse

#note throws unused errors, but it still works...

####Comparison between different sources graph#################################
def CompPlot(comp,ptitle):
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    ax = comp[['SE_gross(kWh/m^3)','SE_net(kWh/m^3)','MW_gross(MW)','MW_net(MW)']] \
        .plot(kind='barh', title =("Source Comparison for Tampa Bay - %s" %ptitle), figsize=(15, 10), legend=True, fontsize=12)
        #plt.legend(['SE_gross(kWh/m^3)','SE_net(kWh/m^3)','MW_gross(MW)','MW_net(MW)'])
    return plt.show()

####Stacked bar graph for Specific Energy######################################
def StackPlotSE(comp,ptitle):
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    N = 3
    E1 = comp['SE_net(kWh/m^3)']
    E2 = comp['SE_Ineff(kWh/m^3)'] + E1
    E3 = comp['SE_PT(kWh/m^3)'] + E2
    E4 = comp['SE_Tr(kWh/m^3)'] +E3
    
    ind = np.arange(N) + .15 # the x locations for the groups
    width = 0.35       # the width of the bars
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, E4, width, color='g') 
    rects2 = ax.bar(ind, E3, width, color='r') 
    rects3 = ax.bar(ind, E2, width, color='b')
    rects4 = ax.bar(ind, E1, width, color='cyan')
    
    # add some text for labels, title and axes ticks
    ax.set_ylabel('Specific Energy (kWh m^-3 mixed solution)')
    #ax.set_title('Source Comparison for Tampa Bay - %s' %ptitle)
    
    ax.set_xticks(ind)
    ax.set_xticklabels( ( 'SW-WW', 'ROC-WW', 'ROC-SW') )
    plt.legend(['Transmission Energy Use','Pretreatment Energy Use','Energy Loss (Inefficiency)','Net Energy Generation'])
    return plt.show()

####Stacked bar graph for Cost#################################################
def StackPlotCost(comp,ptitle):
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    N = 3

    E1 = comp['Cost_Memb($)']
    E2 = comp['Cost_Tr($)'] + E1
    E3 = comp['Cost_Turb($)'] + E2
    E4 = comp['Cost_PT($)'] + E3
    E5 = comp['Cost_Cons($)']
    E6 = comp['Cost($)']

    ind = np.arange(N) + .15 # the x locations for the groups
    width = 0.35       # the width of the bars
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, E6/1e6, width, color='orange') 
    rects2 = ax.bar(ind, E5/1e6, width, color='indigo') 
    rects3 = ax.bar(ind, E4/1e6, width, color='r') 
    rects4 = ax.bar(ind, E3/1e6, width, color='b') 
    rects5 = ax.bar(ind, E2/1e6, width, color='g')
    rects6 = ax.bar(ind, E1/1e6, width, color='cyan')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Cost($) (Millions)')
    #ax.set_title('Cost Comparison for Tampa Bay - %s' %ptitle)
    
    ax.set_xticks(ind)
    ax.set_xticklabels( ( 'SW-WW', 'ROC-WW', 'ROC-SW') )
    plt.legend(['Eng-Legal-Admin Cost','Pipe-Elec-Prep Cost','Pretreatment Cost','Turbine Cost','Transmission Cost','Membrane Cost'])
    return plt.show()


####Stacked bar graph for ROI#################################################
def StackPlotROI(comp,ptitle):
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    N = 3

    E1 = comp['Cost($)']
    E2 = comp['PV_rev($)']

    
    ind = np.arange(N) + .15 # the x locations for the groups
    width = 0.35       # the width of the bars
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, E2/1e6, width, color='g')
    rects2 = ax.bar(ind+width, E1/1e6, width, color='r')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('$ (Millions)')
    #ax.set_title('Potential ROI Comparison for Tampa Bay - %s' %ptitle)
    
    ax.set_xticks(ind)
    ax.set_xticklabels( ( 'SW-WW', 'ROC-WW', 'ROC-SW') )
    plt.legend(['Net Present Value', 'Total Cost'])
    return plt.show()

####Unit Stacked bar graph for Cost############################################
def StackPlotCostUnit(comp,ptitle):
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    N = 3

    E1 = comp['Cost_Memb($)']/(comp['MW_net(MW)'] * 1000)
    E2 = (comp['Cost_Tr($)'] )/(comp['MW_net(MW)'] * 1000) + E1
    E3 = (comp['Cost_Turb($)'] )/(comp['MW_net(MW)'] * 1000) + E2
    E4 = (comp['Cost_PT($)'] )/(comp['MW_net(MW)'] * 1000) + E3
    E5 = comp['Cost_Cons($)']/(comp['MW_net(MW)'] * 1000)
    E6 = comp['Cost($)']/(comp['MW_net(MW)'] * 1000)

    ind = np.arange(N) + .15 # the x locations for the groups
    width = 0.35       # the width of the bars
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, E6/1e3, width, color='orange') 
    rects2 = ax.bar(ind, E5/1e3, width, color='indigo') 
    rects3 = ax.bar(ind, E4/1e3, width, color='r') 
    rects4 = ax.bar(ind, E3/1e3, width, color='b') 
    rects5 = ax.bar(ind, E2/1e3, width, color='g')
    rects6 = ax.bar(ind, E1/1e3, width, color='cyan')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Unit Cost ($ kW^-1) (Thousands)')
    #ax.set_title('Cost Comparison for Tampa Bay - %s' %ptitle)
    
    ax.set_xticks(ind)
    ax.set_xticklabels( ( 'SW-WW', 'ROC-WW', 'ROC-SW') )
    plt.legend(['Eng-Legal-Admin Cost','Pipe-Elec-Prep Cost','Pretreatment Cost','Turbine Cost','Transmission Cost','Membrane Cost'])
    return plt.show()


####Unit Stacked bar graph for ROI#############################################
def StackPlotROIUnit(comp,ptitle):
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    N = 3

    E1 = comp['Cost($)']/(comp['MW_net(MW)'] * 1000)
    E2 = comp['PV_rev($)']/(comp['MW_net(MW)'] * 1000)

    
    ind = np.arange(N) + .15 # the x locations for the groups
    width = 0.35       # the width of the bars
    
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, E2/1e3, width, color='g')
    rects2 = ax.bar(ind+width, E1/1e3, width, color='r')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('$ kW^-1 (Thousands)')
    #ax.set_title('Potential ROI Comparison for Tampa Bay - %s' %ptitle)
    
    ax.set_xticks(ind)
    ax.set_xticklabels( ( 'SW-WW', 'ROC-WW', 'ROC-SW') )
    plt.legend(['Unit Revenue', 'Unit Cost'])
    return plt.show()
####Cost vs Net Power Scatter Plot#############################################
def ScatmPlot(c1,c2,c3,xx,yy,yz,titl,Spec,Spec_inc,line,SysCon2_inc):
    #Spec is a set of specific points to focus on
    #Spec inc is whether or not to include those points
    #line variable indicates whether or not to have a line
    #SysCon2_inc variable indicates whether or not to include trials where SysCon2 = 0. If 
    #SysCon2_inc = 1, then they will be included, if not, then they will be filtered out.
    if SysCon2_inc == 0:
        c1 = c1[c1.SysCon2 == 1]
        c2 = c2[c2.SysCon2 == 1]
        c3 = c3[c3.SysCon2 == 1]
    
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    fig = plt.figure()
    
    #SW vs WW
    x = c1[xx]
    y = c1[yy]/(c1[yz])
    # Fit with polyfit
    b, m = polyfit(x, y, 1)
    plt.plot(x, y, '.', color = 'g')
    if line == 1:
        plt.plot(x, b + m * x, '-', color = 'g')
    
    #ROC vs WW
    x = c2[xx]
    y = c2[yy]/(c2[yz])
    # Fit with polyfit
    b, m = polyfit(x, y, 1)
    plt.plot(x, y, '.', color = 'r')
    if line == 1:
        plt.plot(x, b + m * x, '-', color = 'r')
        
    #ROC vs SW
    x = c3[xx]
    y = c3[yy]/(c3[yz])
    # Fit with polyfit
    b, m = polyfit(x, y, 1)
    plt.plot(x, y, '.', color = 'b')
    if line == 1:
        plt.plot(x, b + m * x, '-', color = 'b')
        
    if Spec_inc == 1:
        x1 = Spec[xx][0]
        x2 = Spec[xx][1]
        x3 = Spec[xx][2]
        y1 = Spec[yy][0]/Spec[yz][0]
        y2 = Spec[yy][1]/Spec[yz][1]
        y3 = Spec[yy][2]/Spec[yz][2]
        plt.scatter(x1, y1, marker = (5,2), s = 40, color = 'g')
        plt.scatter(x2, y2, marker = (5,2), s = 40, color = 'r')
        plt.scatter(x3, y3, marker = (5,2), s = 40, color = 'b')
        
    
    #Chart Details
    #plt.title(titl)
    plt.xlabel("Net SE (kWh m^-3 of mixed solution)")
    plt.ylabel("Unit NPV ($ L^-1 hr of mixed solution)")
    if line == 1:
        plt.legend(['SW-WW','SW-WW','ROC-WW','ROC-WW', 'ROC-SW', 'ROC-SW'])
    else:
        plt.legend(['SW-WW','ROC-WW','ROC-SW'])
    return plt.show()

###Return on Investment over time plot#########################################

def ROITPlot(proj,comp,case,OpTime):
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    #comp is the dataset, case is the gradient of interest (e.g. sw_ww)
    MW_net = comp.loc[case]["MW_net(MW)"]
    Cost = comp.loc[case]["Cost($)"]
    memb = comp.loc[case]["Membrane"]
    Cost_OM = comp.loc[case]["Cost_OM($/Yr)"]
    Proj_vec = []
    Year_vec = []
    Cost_vec = []
    N = proj[0]
    
    for i in range(0,N):
        fpg = pse.finproj(i,proj[1],proj[2],proj[3],MW_net,OpTime,Cost,Cost_OM)
        Proj_vec.append((fpg[0]/1e6))
        Year_vec.append(i)
        Cost_vec.append(Cost/1e6)
        if Proj_vec[i] >= (Cost/1e6):
            break
    
    N = len(Proj_vec)
    ind = np.arange(N)    # the x locations for the groups
    width = 0.35       # the width of the bars: can also be len(x) sequence
    
    fig = plt.figure()
    p1 = plt.bar(ind, Cost_vec, width)
    p2 = plt.bar(ind, Proj_vec, width)
    
    plt.ylabel('$ (Millions)')
    plt.xlabel('Year')
    titl = ('ROI for %s, using the %s Membrane'%(case,memb))
    #plt.title(titl)
    plt.legend(('Capital Cost', 'Return'))
    
    return plt.show()

"""
import pandas as pd
import matplotlib.pyplot as plt
from pandas.plotting import parallel_coordinates


df = pd.read_csv('SA tests/Analyze/Combo.csv')
# Make the plot
fig = plt.figure()
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
parallel_coordinates(df, 'Interval', colormap=plt.get_cmap("viridis"))
plt.grid(b=None)
plt.show()
"""
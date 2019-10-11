# -*- coding: utf-8 -*-
"""
Created on Tue Aug 21 11:04:34 2018

membrane data storage file

@author: Drizdar
"""
import props
import pandas as pd

#%%#Membrane Types#############################################################
# order is: name,id, k [L m-2 h-1], D [m^2/h], S [m], B [L m-2 h-1], A [L m-2 h-1 bar-1], material, ref
#note: standard membrane area of 372 m^2 is assumed - this is basically a legacy product now, but will be kept just in case
#note: temperatures standardized at 25*C

membs =[
    #note: FO flat sheet (cartridge) membrane, parameters defined in Wang, 2010
    #membrane is believed to be made of cellulose acetate
    props.MembraneType("G_06_1","Gray, 2006 - HTI",10.08,5.796e-6,5.75e-4,0.46,1.13,"CTA","Gray, G., McCutcheon, J., Elimelech, M., 2006. Internal concentration polarization in forward osmosis: role of membrane orientation. Desalination 197, 1–8.")
    #note: Membranes from Achilli, 2009 use these standards: B=0.40, A=0.673
    ,props.MembraneType("A_09_1","Achilli, 2009 - #1",303.5,5.299e-6,7.07e-4,0.40,0.673,"CTA","Achilli, A., Cath, T.Y., Childress, A.E., 2009. Power generation with pressure retarded osmosis: An experimental and theoretical investigation. Journal of membrane science.")
    ,props.MembraneType("A_09_2","Achilli, 2009 - #2",305.3,5.364e-6,6.82e-4,0.40,0.673,"CTA","Achilli, A., Cath, T.Y., Childress, A.E., 2009. Power generation with pressure retarded osmosis: An experimental and theoretical investigation. Journal of membrane science.")
    ,props.MembraneType("A_09_3","Achilli, 2009 - #3",311.8,5.602e-6,6.46e-4,0.40,0.673,"CTA","Achilli, A., Cath, T.Y., Childress, A.E., 2009. Power generation with pressure retarded osmosis: An experimental and theoretical investigation. Journal of membrane science.")
    #note: Km calculated as D/S, B and A have +/- ranges, designed for FO
    ,props.MembraneType("C_10_1","Chou, 2010 - FO-#C-TFC",10.54,5.796e-6,5.50e-4,0.22,3.5,"TFC","Chou, S., Shi, L., Wang, R., Tang, C., Qiu, C., Fane, A., 2010. Characteristics and potential applications of a novel forward osmosis hollow fiber membrane. Desalination 261, 365–372.")                          
    #note standard D of 5.796e-6, membranes are hollow fiber, designed for FO
    ,props.MembraneType("W_10_1","Wang, 2010 - FO-#A-TFC",3.96,5.796e-6,1.37e-3,0.29,0.94,"TFC","Wang, R., Shi, L., Tang, C., Chou, S., Qiu, C., Fane, A., 2010. Characterization of novel forward osmosis hollow fiber membranes. J Membrane Sci 355, 158–167.")                          
    ,props.MembraneType("W_10_2","Wang, 2010 - FO-#B-TFC",9.00,5.796e-6,5.95e-4,0.20,2.23,"TFC","Wang, R., Shi, L., Tang, C., Chou, S., Qiu, C., Fane, A., 2010. Characterization of novel forward osmosis hollow fiber membranes. J Membrane Sci 355, 158–167.")                          
    #note membranes from Yip, 2011 use a standard D of 5.364e-6, and a standard k of 138.6+-73.9. 
    #138.6 used for all k values, since individual k values not given.
    ,props.MembraneType("Y_11_1","Yip, 2011 - LP#1",138.6,5.364e-6,3.07e-4,0.16,1.74,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")
    ,props.MembraneType("Y_11_2","Yip, 2011 - LP#2",138.6,5.364e-6,3.55e-4,0.08,1.42,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")
    ,props.MembraneType("Y_11_3","Yip, 2011 - LP#3",138.6,5.364e-6,3.84e-4,0.09,1.71,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")
    ,props.MembraneType("Y_11_4","Yip, 2011 - MP#1",138.6,5.364e-6,3.70e-4,0.88,5.81,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")
    ,props.MembraneType("Y_11_5","Yip, 2011 - MP#2",138.6,5.364e-6,3.32e-4,0.77,4.08,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")    
    ,props.MembraneType("Y_11_6","Yip, 2011 - MP#3",138.6,5.364e-6,3.16e-4,0.61,3.16,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")
    ,props.MembraneType("Y_11_7","Yip, 2011 - HP#1",138.6,5.364e-6,3.27e-4,5.45,7.55,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")
    ,props.MembraneType("Y_11_8","Yip, 2011 - HP#2",138.6,5.364e-6,3.36e-4,4.12,7.35,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")
    ,props.MembraneType("Y_11_9","Yip, 2011 - HP#3",138.6,5.364e-6,4.16e-4,3.86,7.76,"TFC","Yip, N.Y., Tiraferri, A., Phillip, W.A., Schiffman, J.D., Hoover, L.A., Kim, Y.C., Elimelech, M., 2011. Thin-film composite pressure retarded osmosis membranes for sustainable power generation from salinity gradients. Environ. Sci. Technol. 45, 4360–9.")                          
    #note: Diffusion constant copied from Straub, 2014, membrane is hollow fiber
    ,props.MembraneType("C_12_1","Chou, 2012 - PRO-TFC-PES",12.2,5.328e-6,4.60e-4,0.14,3.32,"TFC","Chou, S., Wang, R., Shi, L., She, Q., Tang, C., Fane, A., 2012. Thin-film composite hollow fiber membranes for pressure retarded osmosis (PRO) process with high power density. J Membrane Sci 389, 25–33.")                          
    #note, membranes from She, 2012 have +/- ranges for all paramters
    #there was a 4th membrane listed in the paper, but it was at 35C, not 25C
    ,props.MembraneType("S_12_1","She, 2012 - CTA-P",12.06,5.796e-6,4.80e-4,0.63,0.749,"CTA","She, Q., Jin, X., Tang, C., 2012. Osmotic power production from salinity gradient resource by pressure retarded osmosis: Effects of operating conditions and reverse solute diffusion. J Membrane Sci 401, 262–273.")
    ,props.MembraneType("S_12_2","She, 2012 - CTA-NW",4.21,5.796e-6,1.38e-3,0.07,0.436,"CTA","She, Q., Jin, X., Tang, C., 2012. Osmotic power production from salinity gradient resource by pressure retarded osmosis: Effects of operating conditions and reverse solute diffusion. J Membrane Sci 401, 262–273.")
    ,props.MembraneType("S_12_3","She, 2012 - CTA-W",9.83,5.796e-6,5.90e-4,0.28,0.367,"CTA","She, Q., Jin, X., Tang, C., 2012. Osmotic power production from salinity gradient resource by pressure retarded osmosis: Effects of operating conditions and reverse solute diffusion. J Membrane Sci 401, 262–273.")
    #note: k calculated as D/S, D copied from Chou, 2012, A,S,B have +/- ranges
    ,props.MembraneType("C_13_1","Chou, 2013 - PRO-TFC-PEI",8.73,5.328e-6,6.1e-4,0.24,1.52,"TFC","Chou, S., Wang, R., Fane, A., 2013. Robust and High performance hollow fiber membranes for energy harvesting from salinity gradients by pressure retarded osmosis. Journal of Membrane Science 448, 44–54.")                           
    #note: mass transfer coeff copied from Song, 2013 - membrane is hollow fiber
    ,props.MembraneType("H_13_1","Han, 2013 - TFC-HF1",76.7,5.210e-6,9.87e-4,0.13,1.40,"TFC","Han, G., Wang, P., Chung, T.-S., 2013. Highly Robust Thin-Film Composite Pressure Retarded Osmosis (PRO) Hollow Fiber Membranes with High Power Densities for Renewable Salinity-Gradient Energy Generation. Environ Sci Technology 47, 8070–8077.")
    ,props.MembraneType("H_13_2","Han, 2013 - TFC-HF2",76.7,5.321e-6,7.45e-4,0.41,1.70,"TFC","Han, G., Wang, P., Chung, T.-S., 2013. Highly Robust Thin-Film Composite Pressure Retarded Osmosis (PRO) Hollow Fiber Membranes with High Power Densities for Renewable Salinity-Gradient Energy Generation. Environ Sci Technology 47, 8070–8077.")
    ,props.MembraneType("H_13_3","Han, 2013 - TFC-HF3",76.7,5.331e-6,7.76e-4,0.48,1.90,"TFC","Han, G., Wang, P., Chung, T.-S., 2013. Highly Robust Thin-Film Composite Pressure Retarded Osmosis (PRO) Hollow Fiber Membranes with High Power Densities for Renewable Salinity-Gradient Energy Generation. Environ Sci Technology 47, 8070–8077.")
    #note: A,S,B all have +/- values, K and D are same for all membranes
    #membrane type is a thin-film nanofiber composite, classified as TFC here
    ,props.MembraneType("S_13_1","Song, 2013 - TNC-1",76.7,5.328e-6,1.49e-4,0.28,1.23,"TFC","Song, X., Liu, Z., Sun, D., 2013. Energy recovery from concentrated seawater brine by thin-film nanofiber composite pressure retarded osmosis membranes with high power density. Energy Environ Sci 6, 1199–1210.")
    ,props.MembraneType("S_13_2","Song, 2013 - TNC-2",76.7,5.328e-6,1.35e-4,1.19,3.82,"TFC","Song, X., Liu, Z., Sun, D., 2013. Energy recovery from concentrated seawater brine by thin-film nanofiber composite pressure retarded osmosis membranes with high power density. Energy Environ Sci 6, 1199–1210.")
    ,props.MembraneType("S_13_3","Song, 2013 - TNC-3",76.7,5.328e-6,1.40e-4,3.86,5.31,"TFC","Song, X., Liu, Z., Sun, D., 2013. Energy recovery from concentrated seawater brine by thin-film nanofiber composite pressure retarded osmosis membranes with high power density. Energy Environ Sci 6, 1199–1210.")
    #note, Membrane from Achilli, 2014 is an Oasys Water 4040 spiral wound membrane module.
    #Actual area is 4.18 m^2. Diffusion and mass transfer coeff copied from Straub, 2014.
    ,props.MembraneType("A_14_1","Achilli, 2014",99,5.328e-6,3.10e-4,0.09,5.11,"TFC","Achilli, A., Prante, J.L., Hancock, N.T., Maxwell, E.B., Childress, A.E., 2014. Experimental results from RO-PRO: a next generation system for low-energy desalination. Environ. Sci. Technol. 48, 6437–43.")
    #everything from Straub is actually from the paper.
    ,props.MembraneType("S_14_1","Straub, 2014", 99,5.328e-6,5.64e-4,0.39,2.49,"TFC","Straub, A., Yip, N., Elimelech, M., 2014. Raising the Bar: Increased Hydraulic Pressure Allows Unprecedented High Power Densities in Pressure-Retarded Osmosis. Environ Sci Technology Lett 1, 55–59.")
    #note, membrane in Altaee, 2016 was based on membranes from Achilli, 2009. 
    #Therefore, CTA is assumed to be the membrane material. 
    #Diffusion value is average of 6.1-6.4E-6
    ,props.MembraneType("A_16_1","Altaee, 2016",306,6.250e-6,8.00e-4,0.40,0.67,"CTA","Altaee, A., Millar, G., Zaragoza, G., 2016. Integration and optimization of pressure retarded osmosis with reverse osmosis for power generation and high efficiency desalination. Energy 103, 110–118.")
    ]

#%%#Pump Data##################################################################

Pump_dat = {
    #values are ID, name, W_n, a, n
    #cost curves from "McGivney, W., Kawamura, S., 2008. Cost Estimating Manual for Water Treatment Facilities. John Wiley & Sons, Inc."
    #cost curves cover entire pump station cost, based on $/horsepower
    #cost estimates are class 5 (+50%, -30%) - based on ENR CCI = 8889, Los Angeles California, (April, 2007)
    #Pump Efficiency values from "Hammer, M.J., Hammer, Jr., M.J., 2012. Water and Wastewater Technology, 7th ed. Pearson Education Inc., Upper Saddle River, NJ."
    "Sub":props.PumpType("Sub","Submersible", 0.89, 16599, 0.695) #eff based on river water pump from loeb 2002 paper
    ,"Horiz":props.PumpType("Horiz","Horizontal", 0.75, 62040, 0.6128)
    ,"Vert":props.PumpType("Vert","Submersible", 0.75, 26238, 0.619)
    }

#%%#Pretreatment Data##########################################################

Pt_dat = {
    #note: format is ID,name,PT_E(kWh/m^3),ceqtype, a, b, c, n, X_unit, conv)
    #ceqtype can either be Linear (aX + b), Power (aX^n), or Polynomial (aX^2 + bX + c)
    #if zeros, data unknown or not necessary at the moment
    #energy data from "EPRI, 2013. Electricity Use and Management in the Municipal Water Supply and Wastewater Industries ( No. 3002001433). Electric Power Research Institute."
    #original figures from EPRI reference were in kWh/Mgal, on pg 73
    #unless stated otherwise, cost data from "McGivney, W., Kawamura, S., 2008. Cost Estimating Manual for Water Treatment Facilities. John Wiley & Sons, Inc."
    #cost estimates are class 5 (+50%, -30%) - based on ENR CCI = 8889, Los Angeles California, (April, 2007)
    #conv is number that X is multiplied by, ex, if X is Q, assume going from L/hr to X_unit
    "RswP":props.PTType("RswP","Raw surface water pumping", 0.038,"Linear",12169,60716,0,0,"MGD",6.3401E-6) #X range: 1 - 200, unit represents pump capacity
    ,"RgwP":props.PTType("RgwP","Raw groundwater pumping", 0.243,"Linear",12169,60716,0,0,"MGD",6.3401E-6) #X range: 1 - 200, unit represents pump capacity
    ,"RM":props.PTType("RM","Rapid mixing", 0.011,"N/A",0,0,0,0,"Gal",0) #Unit in book is for basin Vol
    ,"Floc":props.PTType("Floc","Flocculation", 0.003,"N/A",0,0,0,0,"MG",0) #Unit in book is for basin Vol
    ,"Sed":props.PTType("Sed","Sedimentation", 0.004,"N/A",0,0,0,0,"N/A",0)
    ,"Chfs":props.PTType("Chfs","Chemical feed systems", 0.017,"N/A",0,0,0,0,"N/A",0)
    #cost values for MF and UF from "AWWA, 2008. Microfiltration and ultrafiltration membranes for drinking water. Journal (American Water Works Association) 100, 84–97."
    #estimates converted from December 2004 values to September 2007 values using inflation calculator at https://www.bls.gov/data/inflation_calculator.htm
    #lower range of estimate used for MF, upper range for UF
    #note: alternate source used for MF/UF due to lack of individual MF/UF curves in McGivney Source
    ,"MF":props.PTType("MF","Microfiltration (in lieu of sedimentation)", 0.026,"Linear",0.6946,0,0,0,"L/hr",1)
    ,"UF":props.PTType("UF","Ultrafiltration (contaminant removal)", 0.211,"Linear",1.0419,0,0,0,"L/hr",1)
    ,"ROb":props.PTType("ROb","Reverse Osmosis (brackish water)", 1.585,"N/A",0,0,0,0,"N/A",0)
    ,"ROo":props.PTType("ROo","Reverse Osmosis (ocean water)", 3.170,"N/A",0,0,0,0,"N/A",0)
    ,"Daf":props.PTType("Daf","Dissolved air floatation", 0.029,"N/A",0,0,0,0,"N/A",0)
    ,"Ast":props.PTType("Ast","Air stripping", 0.099,"N/A",0,0,0,0,"N/A",0)
    ,"Bwp":props.PTType("Bwp","Backwash water pumps", 0.004,"N/A",0,0,0,0,"SF",0) #Unit in book is for filter surface area
    ,"Rp":props.PTType("Rp","Residuals pumping", 0.001,"N/A",0,0,0,0,"N/A",0)
    ,"DisCl":props.PTType("DisCl","Onsite chlorine for generation", 0.022,"N/A",0,0,0,0,"N/A",0)
    ,"DisOz":props.PTType("DisOz","Ozone disinfection", 0.037,"N/A",0,0,0,0,"N/A",0)
    ,"DisUV":props.PTType("DisUV","UV disinfection", 0.016,"N/A",0,0,0,0,"N/A",0)
    ,"Fwp":props.PTType("Fwp","Finished water pumping", 0.275,"Linear",18888,140743,0,0,"MGD",6.3401E-6) #Note, TDH = 100, X range is from 1.45 to 300, Xunit is for pump capacity
    ,"Npl":props.PTType("Npl","Nonprocess loads (buildings, HVAC, lighting, computers, etc.", 0.079,"Power",63568,0,0,0.553,"MGD",6.3401E-6) #Xrange is 1 - 200 MGD for Plant Capacity, in book, this is for "Administration, Laboratory, and Maintenance Building" 
    }

#%%#Hydraulics Data############################################################
ks_vals = {
        #note: values in m, data from: "Young, D.F., Munson, B.R., Okiishi, T.H., Huebsch, W.W., 2012. Introduction to Fluid Mechanics, 5th ed. John Wiley & Sons, Inc., Hoboken, NJ.",
        #"Chadwick, A., Morfett, J., 1998. Hydraulics in Civil and Environmental Engineering, 3rd ed. E & FN Spon, London.",
        #Engineering ToolBox, (2003). Ventilation Ducts - Roughness & Surface Coefficients. [online] Available at: https://www.engineeringtoolbox.com/surface-roughness-ventilation-ducts-d_209.html [Accessed 22/8/2018]"
        'Riveted steel':(0.9 + 9.0)/2, #Young, 2012
        'Commercial steel':0.045, #Young, 2012
        'Brass':0.003, #Chadwick, 1998
        'Copper':0.003, #Chadwick, 1998
        'Cast iron':0.26, #Young, 2012
        'Wrought iron':0.06, #Chadwick, 1998
        'Galvanized iron':0.15, #Young, 2012
        'Bitumen-lined ductile iron':0.03, #Chadwick, 1998
        'Spun concrete-lined ductile iron':0.03, #Chadwick, 1998
        'Concrete':(0.3 + 3.0)/2, #Young, 2012
        'Slimed concrete sewer':6.0, #chadwick, 1998
        'PVC and Plastic Pipes':(0.0015 + 0.007)/2, #Engineering Toolbox, 2003
        'Drawn tubing': 0.0015 #Young, 2012
        }

kL_vals = {
        #note: values in mm, data from: "Young, D.F., Munson, B.R., Okiishi, T.H., Huebsch, W.W., 2012. Introduction to Fluid Mechanics, 5th ed. John Wiley & Sons, Inc., Hoboken, NJ."
        'Entrance':{
                'Reentrant':0.8,
                'Sharp-edged':0.5,
                'Slightly rounded':0.2,
                'Well rounded':0.04},
        'Exit':{
                'Reentrant':1.0,
                'Sharp-edged':1.0,
                'Slightly rounded':1.0,
                'Well rounded':1.0},
        'Elbows':{
                'Regular 90°, flanged':0.3,
                'Regular 90°, threaded':1.5,
                'Long radius 90°, flanged':0.2,
                'Long radius 90°, threaded':0.7,
                'Long radius 45°, flanged':0.2,
                'Long radius 45°, threaded':0.4},
        '180° return bends':{
                '180° return bend, flanged':0.2,
                '180° return bend, threaded':1.5},
        'Tees':{
                'Line flow, flanged':0.2,
                'Line flow, threaded':0.9,
                'Branch flow, flanged':1.0,
                'Branch flow, threaded':2.0},
        'Union,threaded':0.08,
        'Valves':{
                'Globe, fully open':10,
                'Angle, fully open':2,
                'Gate, fully open':0.15,
                'Gate, 1/4 closed':0.26,
                'Gate, 1/2 closed':2.1,
                'Gate, 3/4 closed':17,
                'Swing check, forward flow': 2,
                'Swing check, backward flow':1e9, #infinity
                'Ball valve, fully open':0.05,
                'Ball valve, 1/3 closed':5.5,
                'Ball valve, 2/3 closed':210}
        }

C_Vals = {
        #note: C vals for Hazen-WIlliams Equation, data from: "Hammer, M.J., Hammer, Jr., M.J., 2012. Water and Wastewater Technology, 7th ed. Pearson Education Inc., Upper Saddle River, NJ."
        'New riveted steel':110, #Hammer, 2012
        'New welded steel':120, #Hammer, 2012
        'Copper':(130 + 140)/2, #Hammer, 2012
        '20 year old unlined ductile iron':100, #Hammer, 2012
        '5 year old unlined ductile iron':120, #Hammer, 2012
        'New unlined ductile iron':130, #Hammer, 2012
        'Cement-lined ductile iron':(130 + 150)/2, #Hammer, 2012
        'Concrete':130, #Hammer, 2012
        'Plastic':(140 + 150)/2, #Hammer, 2012   
        }

#%%#Spacer Data############################################################
Spacer_dat = {
    #values are shape,alpha,beta,a_hc,a_b,shape_type (0 = circle, 1 = ellipse, 2 = square, 3 = wing)
    #values from "Guillen, G., Hoek, E., 2009. Modeling the impacts of feed spacer geometry on reverse osmosis and nanofiltration processes. Chemical Engineering Journal."
    #Table 3 in aforementioned reference
    "1:3 ellipse":props.SpacerType("1:3 ellipse", 0.13, 128.24, 0.13, (1/3), 1)
    ,"1:2 ellipse":props.SpacerType("1:2 ellipse", 0.17, 141.24, 0.17, 0.5, 1)
    ,"1:3 wing":props.SpacerType("1:3 wing", 0.19, 135.56, 0.18, (1/3), 3)
    ,"Square":props.SpacerType("Square", 0.42, 176.99, 0.20, 1, 2)
    ,"1:2 wing":props.SpacerType("1:2 wing", 0.29, 152.16, 0.22, 0.5, 3)
    ,"Circle":props.SpacerType("Circle", 0.42, 189.29, 0.25, 1, 0)
    ,"1:1 wing":props.SpacerType("1:1 wing", 0.59, 190.96, 0.27, 1, 3)
    ,"2:1 wing":props.SpacerType("2:1 wing", 0.91, 236.73, 0.30, 2, 3)
    ,"2:1 ellipse":props.SpacerType("2:1 ellipse", 1.44, 315.63, 0.33, 2, 1)
    }
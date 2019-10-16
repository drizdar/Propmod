# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 15:06:17 2018

This file contains classes for the various plants and devices that are used 
for PRO. Really it's just a holding file for classes haha.

@author: Drizdar
"""
import pandas as pd
import math

#%%#Plants#####################################################################
class DesalPlant:
   'Common base class for different desalination plants'

   def __init__(self, name, C_ROC, C_SW, Q_ROC, OpTime, Elev):
      self.name = name
      self.C_ROC = C_ROC #g/kg
      self.C_SW = C_SW #g/kg
      self.Q_ROC = Q_ROC #m^3/s
      self.OpTime = OpTime #% of time plant operates (in decimal form)
      self.thryr = OpTime * 24 * 365 
      self.Elev = Elev #m above sea level
   
   def displayStats(self):
      print("Name: ", self.name,  
            ", Brine Concentration:",self.C_ROC,"g/kg",
            ", Saltwater Concentration:",self.C_SW,"g/kg",
            ", Maximum Brine Flow Rate:",self.Q_ROC,"m^3/s",
            ", Average Yearly Operating Time:",self.thryr,"hrs per year",
            ", Elevation:",self.Elev,"m above sea level")
            


class WWPlant:
   'Common base class for different wastewater plants'

   def __init__(self, name, C_Eff, C_SW, Q_Eff, OpTime, Elev):
      self.name = name
      self.C_Eff = C_Eff #g/kg
      self.C_SW = C_SW #g/kg
      self.Q_Eff = Q_Eff #m^3/s
      self.OpTime = OpTime #% of time plant operates (in decimal form)
      self.thryr = OpTime * 24 * 365 
      self.Elev = Elev #m above sea level
   
   def displayStats(self):
      print("Name: ", self.name,  
            ", Brine Concentration:",self.C_Eff,"g/kg",
            ", Saltwater Concentration:",self.C_SW,"g/kg",
            ", Maximum Brine Flow Rate:",self.Q_Eff,"m^3/s",
            ", Average Yearly Operating Time:",self.thryr,"hrs per year",
            ", Elevation:",self.Elev,"m above sea level")

#%%#Devices####################################################################
class TurbineType:
   'Common base class for different turbines'

   def __init__(self, name,W_n):
      self.name = name
      self.W_n = W_n
   
   def displayStats(self):
      print("Name: ", self.name,  
            ", Generation Efficiency:",self.W_n,"%%")

class PumpType:
   'Common base class for different pumps'
   #under construction

   def __init__(self, ID, name, W_n, a, n):
      self.ID = ID
      self.name = name
      self.W_n = W_n
      self.a = a
      self.n = n
   
   def to_dict(self):
      return {
               'ID': self.ID,
               'name': self.name,
               'W_n': self.W_n,
               'a': self.a,
               'n': self.n,
               }
       
   def displayStats(self):
      print("Name: ", self.name,  
            ", Pumping Efficiency:",self.W_n,"%%")
      
class PTType:
   'Common base class for different pretreatment methods'
   'note: pretreatment will be used for SW, ROC and WW are assumed to be clean'

   def __init__(self, ID, name, PT_E, ceqtype, a, b, c, n, X_unit, conv):
      self.ID = ID
      self.name = name
      self.PT_E = PT_E #kWh/m^3
      self.ceqtype = ceqtype
      self.a = a
      self.b = b
      self.c = c
      self.n = n
      self.X_unit = X_unit
      self.conv = conv
      
   def to_dict(self):
       return {
               'ID': self.ID,
               'name': self.name,
               'PT_E': self.PT_E,
               'ceqtype': self.ceqtype,
               'a': self.a,
               'b': self.b,
               'c': self.c,
               'n': self.n,
               'X_unit': self.X_unit
               }
   
   def costQ(self,X,inf):
      X = X * self.conv #converts X from L/hr to X_unit
      if self.ceqtype == 'Polynomial':
          cost = (self.a * X**2 + self.b * X + self.c) * inf
      elif self.ceqtype == 'Power':
          cost = (self.a * X**self.n) * inf
      elif self.ceqtype == 'Linear':
          cost = (self.a * X + self.b) * inf
      else:
          cost = 0
      return cost 
   
   def displayStats(self):
      print("Name: ", self.name,  
            ", Specific Energy:",self.W_n,"kWh/m^3")
      
      
class MembraneType:
   '''
   Common base class for different membrane types.
   '''

   def __init__(self, name,id, k, D, S, B, A, mat, ref):
      self.id = id
      self.name = name
      self.k = k
      self.D = D
      self.S = S
      self.B = B
      self.A = A
      self.mat = mat
      self.ref = ref

   def to_dict(self):
       return {
               '0_Name': self.name,
               'k (L m-2 h-1)': self.k,
               'D (m^2/h)': self.D,
               'S (m)': self.S,
               'B (L m-2 h-1)': self.B,
               'A (L m-2 h-1 bar-1)': self.A,
               '0_mat': self.mat,
               }

   def displayStats(self):
      print("Name:", self.name,  
            ", Mass Transfer Coefficient:",self.k,"L m-2 h-1",
            ", Diffusion Coefficient:",self.D,"m^2/s",
            ", Structural Parameter:",self.S,"m",
            ", Salt Permeability of the Active Layer:",self.B,"L m-2 h-1",
            ", Water Permeability of the Active Layer:",self.A,"L m-2 h-1 bar-1",
            ", Material:",self.mat,
            ", Reference:",self.ref
            )
      
#%%#Options####################################################################
##Put here the potential gradient combination parameters
##e.g., distance, pt use, cn_eff      


class Tr_use:
   'Common base class for transmission energy usage between plants'

   def __init__(self, DesalPlant,WWPlant,ks,rho,D,mu,z1,z2,L,g,Q,Pump,inf,kLtot):
      self.DesalPlant = DesalPlant
      self.WWPlant = WWPlant
      self.ks = ks
      self.rho = rho #kg/m^3
      self.D = D #m
      self.mu =  mu #kg m-1 s-1
      self.z1 = z1 #m above sea level
      self.z2 = z2 #m above sea level
      self.L = L #m
      self.g = g #m/s^2
      self.Q = Q #L/hr
      self.Pump = Pump #%
      self.inf = inf #inflation rate from 2007
      self.kLtot = kLtot
      
   def calc(self):
      import pro_sys_eqns as pse
      Tr_use = pse.engeq(self.ks,self.D,self.rho,self.mu,self.z1,self.z2, \
                         self.L,self.g,self.Q,self.Pump,self.inf,self.kLtot) #kWh/m^3
      Tr_use_dat = pd.Series([Tr_use[0],Tr_use[1],Tr_use[2],Tr_use[3],Tr_use[4],Tr_use[5],Tr_use[6],Tr_use[7],Tr_use[8],Tr_use[9]],
                             index=['Tr_use(kWh/m^3)', 'Pump Head(m)', 'P(kW)', 'D(m)','V(m/s)','Pipe Cost($)','Pump Cost($)','Total Cost($)','a','n'])
      #s1 = pd.Series([1,2,3,4,5,6], index=pd.date_range('20130102', periods=6))
      return Tr_use_dat
   
   def displayStats(self):
      print("The", self.DesalPlant,  
            "and the",self.WWPlant,
            "are",self.L/1000,"kilometers apart.")

##Spacers
class SpacerType:
    'common base class for different spacer geometries'
    
    def __init__(self,shape,alpha,beta,a_hc,a_b,shape_type):
        self.shape = shape
        self.alpha = alpha
        self.beta = beta
        self.a_hc = a_hc #a/h_c #ratio to calculate a based on h_c value
        self.a_b = a_b #a/b - ratio to calculate b based on a value
        self.shape_type = shape_type #(0 = circle, 1 = ellipse, 2 = wing, 3 = square)
    
    def displayStats(self):
        print("Shape:", self.shape,
              ", alpha:", self.alpha,
              ", beta:", self.beta,
              ", effective cross section:", self.a_hc,
              ", a to b ratio:", self.a_b)
    
    def perimeter_calc(self, h_c):
        #calculates filament perimeter
        a = self.a_hc * h_c
        b = a/self.a_b
        if self.shape_type == 0: #circle
            P_f = 2 * math.pi * a
        elif self.shape_type == 1: #ellipse
            P_f = math.pi*(3*(a+b) - math.sqrt((3*a + b)*(a + 3*b)))
        elif self.shape_type == 2: #square
            P_f = 8 * a
        elif self.shape_type == 3: #wing
            P_f = (math.pi * a) + (2 * math.sqrt(a**2 + b**2))
        return P_f
    
    def cx_calc(self, h_c):
        a = self.a_hc * h_c
        b = a/self.a_b
        if self.shape_type == 0: #circle
            X_f = math.pi * (a**2)
        elif self.shape_type == 1: #ellipse
            X_f = math.pi * a * b
        elif self.shape_type == 2: #square
            X_f = (2*a)**2
        elif self.shape_type == 3: #wing
            X_f = (0.5 * math.pi * (a**2)) + ((a * b)/2)
        return X_f
        
    def d_h_calc(self,h_c,l_f):
        P_f = self.perimeter_calc(h_c) #filament perimeter
        X_f = self.cx_calc(h_c) #filament cross-sectional area
        d_h = (4*(l_f * h_c - X_f)) / (2*l_f + P_f)
        return d_h
    
    def to_dict(self,h_c,l_f):
       return {
               'Shape': self.shape,
               'alpha': self.alpha,
               'beta': self.beta,
               'd_h': self.d_h_calc(h_c,l_f)
               }
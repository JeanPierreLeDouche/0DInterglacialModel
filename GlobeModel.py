# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 11:47:29 2020

@author: Gebruiker
"""
import numpy as np

# simulation parameters
dt = 10 # yrs
time = 200 * 1e3 # yrs

# constants
r_e = 6371 * 1e3 # m, earth radius

# coefficients per modelling variable (initial values)
# mass
a_0 = -10**-10
a_1 = 10**-9 
a_2 = 10**16
a_3 = 10**8
a_4 = 10**4
a_5 = 10**8
a_coeffs = np.array(a_0, a_1, a_2, a_3, a_4, a_5)

# insolation
S0 = 1360.8 
lat = 45

# obliquity
eps_min = 0.3847 # rad
eps_max = 0.4277 # rad
eps_period = 41040 # yr
phi0 = -152.4972 # rad

# surface temperature
b_0 = 1
b_1 = 1
b_2 = 1

#ocean temperature 
c_0 = 1

# isostatic depression
d_o = 1

# atmospheric CO2
c_i = 270 # ppm 
d_1 = 1

#other
T_ref = 273 ### ???
T_min = 233

def dmdmt(m, T_s, T_o, D, coef, P_max, T_min, T_r ):
    
    # function that calculates dm/dt, input list of a coefficients as coef
    
    P_sl = P_max - P_max/2 * (T_s - 273.15)
    accum = ( 0.25 + coef[0] * m**(1/3)) * (P_sl + coef[1] * m **(1/3))*coef[2]*m**(2/3)
    surf_abl = -1 * (coef[3] * T_s - coef[4]*m**(1/3))
    mar_abl = -1 * D * (T_o - T_r)**2 
    
    mass_change = accum + surf_abl + mar_abl 
    return mass_change
    
def insolation(t):
    eps = eps_min + (eps_max)*np.sin((2*np.pi*t)/eps_period + phi0)
    S = S0*(np.sin(lat*np.pi/180)*np.sin(eps) + np.cos(lat*np.pi/180)*np.cos(eps))
    return S

def T_surf(m, S_E, CO2, coef):
    T = (coef[0] - coef[1] * m**(2/3))*S_E + coef[2] * np.ln(CO2)
    return T

def dTodt(T_s, T_o, coef):
    temp_change = coef * (T_s - T_o) 
    return temp_change

def dIdt(m, D, coef):
    isost_change = coef* (m/3. - D)
    return isost_change 

def dCO2dt(CO2, T_o, coef):
    # introducing a "realistic minimum temperature of the earth" from the long term record
    T_min2 = 271 # K
    CO2_change = coef * (0.0423/(2000/dt))*CO2*(T_o - T_min2)
    return CO2_change

# initial values 
T_o = 278 # K 
T_s = 280 # K 
m = 0.4 # GMSLR 
D = 1 # m 
I = 240 # should be w/m^2 
CO2 = 330  # ppm

### main

T_o_arr = np.zeros((1, int(time/dt)))
T_s_arr = np.zeros((1, int(time/dt)))
m_arr = np.zeros((1, int(time/dt)))

D_arr = np.zeros((1, int(time/dt)))
I_arr = np.zeros((1, int(time/dt)))
CO2_arr = np.zeros((1, int(time/dt)))

for t in np.arange(0, time, dt): 
                    
    















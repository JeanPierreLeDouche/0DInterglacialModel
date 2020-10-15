# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 11:47:29 2020

@author: Gebruiker
"""
import numpy as np
import matplotlib.pyplot as pl

# simulation parameters
dt = 10 # yrs
time = int(200 * 1e3) # yrs

# constants
r_e = 6371 * 1e3 # m, earth radius

# coefficients per modelling variable (initial values)
# mass
a_0 = -10**-6
a_1 = 10**-5
a_2 = 10**1
a_3 = 10**-1
a_4 = 10**0
a_5 = 10**0
a_coeffs = np.array((a_0, a_1, a_2, a_3, a_4, a_5))

# insolation
S0 = 1360.8 
lat = 45

# obliquity
eps_min = 0.3847 # rad
eps_max = 0.4277 # rad
eps_period = 41040 # yr
phi0 = -152.4972 # rad

# surface temperature
b_0 = 1.0*10**(-1)
b_1 = 10.0**(-8)
b_2 = 2*(10.0**1)
b_coeffs = np.array((b_0, b_1, b_2))
T_s_min = 200.
T_s_max = 300.

#ocean temperature 
c_0 = 5*10.0**(-5)

# isostatic depression
d_0 = 10.0**(-15)

# atmospheric CO2
e_0 = 270. # ppm 
e_1 = 5*10.0**(-7)
e_coeffs = np.array((e_0, e_1))
CO2_min = 200.

#other
P_max = 10**-3 #kg per km^2 per yr
T_ref = 273. ### ???
T_min = 233.

def dmdt(m, T_s, T_o, D, coef):
    
    # function that calculates dm/dt, input list of a coefficients as coef
    if (T_s > T_min) and (T_s < 273.15):
        P_sl = P_max*(T_s - T_min)/(273.15 - T_min)
    elif (T_s > 273.15) and (T_s < 275.15):
        P_sl = P_max - (P_max*(T_s - 273.15)/2)
    else:
        P_sl = 0
    
    accum = ( 0.25 + coef[0] * m**(1/3)) * (P_sl + coef[1] * m **(1/3))*coef[2]*m**(2/3)
    surf_abl = -1 * (coef[3] * T_s - coef[4]*m**(1/3))
    mar_abl = -coef[5] * D * (T_o - T_ref)**2 
    
    mass_change = accum + surf_abl + mar_abl 
    return mass_change
    
def Insol(t):
    # depends on t in years ? 
    eps = eps_min + (eps_max)*np.sin((2*np.pi*t)/eps_period + phi0)
    S = S0*(np.sin(lat*np.pi/180)*np.sin(eps) + np.cos(lat*np.pi/180)*np.cos(eps))
    return S

def T_surf(m, I, CO2, coef):
    T = (coef[0] - coef[1] * m**(2/3))*I + coef[2] * np.log(CO2)
    # term1 = (coef[0] - coef[1] * m**(2/3))*I
    # term2 = coef[2] * np.log(CO2)
    return T

def dTodt(T_s, T_o, coef):
    temp_change = coef * (T_s - T_o) 
    return temp_change

def dDdt(m, D, coef):
    isost_change = coef* (m/3. - D)
    return isost_change 

def dCO2dt(CO2, T_o, e_coeffs):
    # introducing a "realistic minimum temperature of the earth" from the long term record
    T_min2 = 271 # K
    CO2_change = e_coeffs[1] * (0.0423/(2000/dt))*CO2*(T_o - T_min2)
    return CO2_change

# initial values 
T_o = 278. # K 
T_s = 280. # K 
m = 10**7 # Pg ice
D = 1. # m 
CO2 = 330.  # ppm
D = 1 # m 
I = Insol(0)
CO2 = 330  # ppm

### main

# declare output arrays
T_o_arr = np.zeros((1, int(time/dt+1)))
T_s_arr = np.zeros((1, int(time/dt +1)))
m_arr = np.zeros((1, int(time/dt +1)))

D_arr = np.zeros((1, int(time/dt +1)))
I_arr = np.zeros((1, int(time/dt +1)))
CO2_arr = np.zeros((1, int(time/dt +1)))


for t in range(0, time+1, dt): 
    
    # first calculate new values for all model variables using the old values
    T_o_new = T_o + dTodt(T_s, T_o, c_0)*dt
    I_new = Insol(t)
    D_new = D + dDdt(m, D, d_0)*dt                

    CO2_new = CO2 + dCO2dt(CO2, T_o, e_coeffs )*dt
    if CO2_new < CO2_min:
        CO2_new = CO2_min
    
    m_new = m + dmdt(m, T_s, T_o, D, a_coeffs)*dt
    if m_new < 0:
        m_new = 0

    T_s_new = T_surf(m, I, CO2, b_coeffs )
    if T_s_new > T_s_max:
        T_s_new = T_s_max
    if T_s_new <  T_s_min:
        T_s_new = T_s_min                 


    # then when everything is calculated replace the values by the new ones
    
    T_o = T_o_new
    I = I_new
    D = D_new
    
    CO2 = CO2_new 
    m = m_new
    T_s = T_s_new 
    
    # save these values to their respective array
    T_o_arr[0,t//10] = T_o
    T_s_arr[0,t//10] = T_s
    D_arr[0,t//10] = D
    
    I_arr[0,t//10] = I
    CO2_arr[0,t//10] = CO2
    m_arr[0,t//10] = m
        
        







## PLOTS
##TODO: get time axes to show years ago

t_axis_f = np.arange(0,time+1,dt) #year, counting up from starting point
t_axis_r = np.zeros(len(t_axis_f)) #years ago
for j in range(len(t_axis_r)):
    t_axis_r[j] = t_axis_f[-(j+1)]


# ice mass
pl.plot(t_axis_f, m_arr[0,:], color='cyan')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel('ice mass [kg]', fontsize=14)
#pl.title('')
pl.grid(True)
pl.show()

# surface temp
pl.plot(t_axis_f, T_s_arr[0,:], color='red')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel('surface temp [K]', fontsize=14)
#pl.title('')
pl.grid(True)
pl.show()

# insolation
pl.plot(t_axis_f, I_arr[0,:], color='orange')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel('mean insolation [W m^-2]', fontsize=14)
#pl.title('')
pl.grid(True)
pl.show()

# CO_2 concentration
pl.plot(t_axis_f, CO2_arr[0,:], color='black')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel('CO_2 concentration [ppm]', fontsize=14)
#pl.title('')
pl.grid(True)
pl.show()

# ocean temp
pl.plot(t_axis_f, T_o_arr[0,:], color='blue')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel('ocean temp [K]', fontsize=14)
#pl.title('')
pl.grid(True)
pl.show()

# isostatic depression
pl.plot(t_axis_f, D_arr[0,:], color='grey')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel('surface depression [m]', fontsize=14)
#pl.title('')
pl.grid(True)
pl.show()













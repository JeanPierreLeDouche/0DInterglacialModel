# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 11:47:29 2020

@author: Gebruiker
"""
# TO DO:
# find initial values from paleoclimate data



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
a_3 = 10**2
a_4 = 10**0
a_5 = 10**-2
a_6 = 10**-2
a_coeffs = np.array((a_0, a_1, a_2, a_3, a_4, a_5, a_6))

# insolation
S0 = 1360.8 
lat = 45

# obliquity
eps_min = 0.3847 # rad
eps_max = 0.4277 # rad
eps_period = 41040 # yr
phi0 = -152.4972 # rad

# surface temperature
b_0 = 230 #initial guess, better to calculate at some point
b_1 = 10.0**(-1)
b_2 = 10.0**(-4)
b_3 = 1*(10.0**1)
b_coeffs = np.array((b_0, b_1, b_2, b_3))
T_s_min = 200.
T_s_max = 300.

#ocean temperature 
c_0 = 5*10.0**(-4) #might need faster time scale for ocean surface

# isostatic depression
d_0 = 10.0**(-15)

# atmospheric CO2
e_0 = 270. # ppm 
e_1 = 1*10.0**(-2)
e_coeffs = np.array((e_0, e_1))
CO2_min = 100.
CO2_max = 600.

#other
P_max = 10**-3 #Pg per km^2 per yr
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
    
    accum = ( 0.25 + coef[0] * m**(1/3)) * (P_sl + coef[1] * m **(1/3))*coef[2]*(coef[3] + m**(2/3))

    if T_s > T_ref:
        surf_abl = -1 * (coef[4] * (T_s - T_ref)*m**(2/3) - coef[5]*m)
    elif T_s <= T_ref:
        surf_abl=0

    if T_o > T_ref:
        mar_abl = -coef[6] * D * (T_o - T_ref)**2 * m**(2/3)
    elif T_o <= T_ref:
        mar_abl = coef[6] * D * (T_o - T_ref)**2 * m**(2/3)                                        
                                                    
    mass_change = accum + surf_abl + mar_abl 
    return mass_change, accum, surf_abl, mar_abl 
    
def Insol(t):
    # depends on t in years ? 
    eps = eps_min + (eps_max)*np.sin((2*np.pi*t)/eps_period + phi0)
    S = S0*(np.sin(lat*np.pi/180)*np.sin(eps) + np.cos(lat*np.pi/180)*np.cos(eps))
    return S

def T_surf(m, I, CO2, coef):
    T = coef[0] + coef[1]*(1 - coef[2] * m**(2/3))*I + coef[3] * np.log(CO2) #recalc coef0
    
    term1 = coef[0]
    term2 = coef[1]*(1 - coef[2] * m**(2/3))*I
    term3 = coef[3] * np.log(CO2)
    
    return T, term1, term2, term3

def dTodt(T_s, T_o, coef):
    temp_change = coef * (T_s - T_o) 
    return temp_change

def dDdt(m, D, coef):
    isost_change = coef* (m/3. - D)
    return isost_change 

def f_CO2(CO2, T_o, e_coeffs):
    # introducing a "realistic minimum temperature of the earth" from the long term record
    T_min2 = 271 # K
    CO2_change = e_coeffs[1] * 0.0423*CO2*(T_o - T_min2) #TODO: need shorter equi time scale; should use deep ocean temp, not ocean surf
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

zeros = np.zeros((1, int(time/dt +1)))

# declare output arrays
T_o_arr = np.zeros((1, int(time/dt+1)))
T_s_arr = np.zeros((1, int(time/dt +1)))
m_arr = np.zeros((1, int(time/dt +1)))

D_arr = np.zeros((1, int(time/dt +1)))
I_arr = np.zeros((1, int(time/dt +1)))
CO2_arr = np.zeros((1, int(time/dt +1)))

m_abl_arr = np.zeros((1, int(time/dt +1)))
acc_arr = np.zeros((1, int(time/dt +1)))
s_abl_arr = np.zeros((1, int(time/dt +1)))

comp1 = np.zeros(int(time/dt +1))
comp2 = np.zeros(int(time/dt +1))
comp3 = np.zeros(int(time/dt +1))

for t in range(0, time+1, dt): 
    
    # first calculate new values for all model variables using the old values
    T_o_new = T_o + dTodt(T_s, T_o, c_0)*dt
    I_new = Insol(t)
    D_new = D + dDdt(m, D, d_0)*dt                

    CO2_new = f_CO2(CO2, T_o, e_coeffs )*dt
    if CO2_new < CO2_min:
        CO2_new = CO2_min
    if CO2_new > CO2_max:
        CO2_new = CO2_max 
    
    dm_dt = dmdt(m, T_s, T_o, D, a_coeffs)
    
    m_new = m + dm_dt[0]*dt
    if m_new < 0:
        m_new = 0

    T_surface = T_surf(m, I, CO2, b_coeffs )  
    T_s_new = T_surface[0]
    
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
    
    m_abl_arr[0,t//10] = dm_dt[3]
    s_abl_arr[0,t//10] = dm_dt[2]
    acc_arr[0, t//10] = dm_dt[1]
    
    comp1[t//10] = T_surface[1]
    comp2[t//10] = T_surface[2]
    comp3[t//10] = T_surface[3]
    


## PLOTS
##TODO: get time axes to show years ago
##TODO: plot mass change components
##TODO: plot ice mass in GMLSE

t_axis_f = np.arange(0,time+1,dt) #year, counting up from starting point
t_axis_r = np.zeros(len(t_axis_f)) #years ago
for j in range(len(t_axis_r)):
    t_axis_r[j] = t_axis_f[-(j+1)]


# ice mass
pl.plot(t_axis_f, np.log10(m_arr[0,:]), color='cyan', label = 'ice mass')
pl.plot(t_axis_f, np.log10(np.abs( m_abl_arr[0,:])), color = 'blue', label = 'marine ablation')
pl.plot(t_axis_f, np.log10(np.abs( s_abl_arr[0,:])), color = 'yellow', label = 'surface ablation')
pl.plot(t_axis_f, np.log10(acc_arr[0,:]), color = 'orange', label = 'accumulation')

#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel('ice mass [Gt]', fontsize=14)
pl.legend()
#pl.title('')
pl.grid(True)
pl.show()

# surface & upper ocean temp
pl.plot(t_axis_f, T_s_arr[0,:], color='red', label = 'Surface temperature')
pl.plot(t_axis_f, T_o_arr[0,:], color='blue', label = 'Ocean temperature')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel('temperature [K]', fontsize=14)
#pl.title('')
pl.grid(True)
pl.legend()
pl.show()

# insolation
pl.plot(t_axis_f, I_arr[0,:], color='orange')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel(r'peak insolation [$W m^{-2}$]', fontsize=14)
#pl.title('')
pl.grid(True)
pl.show()

# CO_2 concentration
pl.plot(t_axis_f, CO2_arr[0,:], color='black')
#pl.axis([,,,])  # define axes 
#pl.xticks(N.arange(,,), fontsize=12) 
#pl.yticks(N.arange(,,), fontsize=12) 
pl.xlabel('time [years]', fontsize=14)
pl.ylabel(r'$CO_2$ concentration [ppm]', fontsize=14)
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


# diagnostics
pl.plot(t_axis_f, comp1[:], color='r', label = 'component1')
pl.plot(t_axis_f, comp2[:], color='r', label = 'component2')
pl.plot(t_axis_f, comp3[:], color='r', label = 'component3')
pl.grid(True)
pl.show()










# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 11:15:05 2020

@author: Gebruiker
"""

import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as pl
import time as tim 
from smt.sampling_methods import LHS
import random as r
import numpy as np
import matplotlib.pyplot as pl
import time as tim 
import itertools 

# simulation parameters
dt = 10 # yrs
time = int(200 * 1e3) # yrs

n_runs = 100



# constants
r_e = 6371 * 1e3 # m, earth radius
Gt_to_SLE_conv = 1 /( 361. * 1e3 ) # multiply value in Gt to get SLE in m

# coefficients per modelling variable (initial values)
# mass
a_0 = -10**-7
a_1 = 9*10**-6
a_2 = 7*10**0
a_3 = 5*10**4
a_4 = 4*10**-4
a_5 = 1*10**-6
a_6 = 5*10**-7
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
b_0 = 208 #initial guess, better to calculate at some point
b_1 = 2 * 10.0**(-2)
b_2 = 10.0**(-9)
b_3 = 0.9*(10.0**1)
b_coeffs = np.array((b_0, b_1, b_2, b_3))
T_s_min = 200.
T_s_max = 300.

#ocean temperature 
c_0 = 5*10.0**(-4) #might need faster time scale for ocean surface

# isostatic depression
d_0 = 10.0**(-8)
d_1 = 1.5*10.0**4
d_coeffs = np.array((d_0,d_1))

# atmospheric CO2
e_0 = 270. # ppm 
e_1 = 1*10.0**(1)
e_2 = 100.
CO2_min = e_2
CO2_max = 600.
e_coeffs = np.array((e_0, e_1,e_2))

#other
P_max = 3*10**-3 #Pg per km^2 per yr

# P_max = P_max * Gt_to_SLE_conv 

T_ref = 273.15
T_min = 228.15

# function that calculates dm/dt, input list of a coefficients as coef
def dmdt(m, T_s, T_o, D, coef):
    
    # precipitation, three cases
    if (T_s > T_min) and (T_s <= 273.15):
        P_sl = P_max*(T_s - T_min)/(273.15 - T_min)
    elif (T_s > 273.15) and (T_s < 298.15):
        P_sl = P_max - (P_max*(T_s - 273.15)/25)
    else:
        P_sl = 0
    
    accum = ( 0.25 + coef[0] * m**(1/3)) * (P_sl + coef[1] * m **(1/3))*coef[2]*(coef[3] + m**(2/3))

    if T_s > T_ref:
        surf_abl = -1 * np.max([(coef[4] * (T_s - T_ref)*m**(2/3) - coef[5]*m),0])
    elif T_s <= T_ref:
        surf_abl=0
    else: # rucksichtloss geimplementiert, warning (!!!)
        surf_abl = 0

    if T_o > T_ref:
        mar_abl = -coef[6] * D * (T_o - T_ref)**2 * m**(2/3)
    elif T_o <= T_ref:
        mar_abl = coef[6] * D * (T_ref - T_o)**2 * m**(2/3)                                        
    else: 
        mar_abl = 0 #achtung: rucksichtloss geimplementiert (!!!)
                                               
    mass_change = accum + surf_abl + mar_abl 

    return mass_change, accum, surf_abl, mar_abl 
    
def Insol(t):
    # depends on t in years ? 
    eps = eps_min + (eps_max - eps_min)*np.sin((2*np.pi*t)/eps_period + phi0)
    S = S0*(np.sin(lat*np.pi/180)*np.sin(eps) + np.cos(lat*np.pi/180)*np.cos(eps))
    return S

def T_surf(m, I, CO2, coef):
    T = coef[0] + coef[1]*( np.max([1 - coef[2] * m**(2/3),0.2]) )*I + coef[3] * np.log(CO2) #
    
    term1 = coef[0]
    term2 = coef[1]*( np.max([1 - coef[2] * m**(2/3),0.2]) )*I
    term3 = coef[3] * np.log(CO2)
    
    return T, term1, term2, term3

def dTodt(T_s, T_o, coef):
    temp_change = coef * (T_s - T_o) 
    return temp_change

def dDdt(m, D, coef):
    isost_change = coef[0]* (m/3. - coef[1]*D)
    return isost_change 

def f_CO2(CO2, T_o, e_coeffs):
    # introducing a "realistic minimum temperature of the earth" from the long term record
    T_min2 = 271 # K
    CO2_change = e_coeffs[1] *(T_o - T_min2) + e_coeffs[2]
    return CO2_change

def scorefunction(arrays): 
    # assuming arrays is a tuple object with arrays for CO_2, T_s and  ice mass
    
    # exlude spinup time of 50 ka: 50*1e3 /10 = 5000 timesteps
    spinup = 5000
    
    CO2_max = np.max(arrays[0, spinup:])
    CO2_min = np.min(arrays[0, spinup:])
    
    T_s_max = np.max(arrays[1, spinup:])
    T_s_min = np.min(arrays[1, spinup:])
    
    m_max = np.max(arrays[2, spinup:]*Gt_to_SLE_conv)
    m_min = np.min(arrays[2, spinup:]*Gt_to_SLE_conv)
    
    # write a bunch of if statements here 
    points = 0
    
    if CO2_max < 300 and CO2_min> 180:
        points += 5
        print('full score for CO2 !')
    elif 300 < CO2_max < 350 and 150 < CO2_min < 180:
        points += 2 
    elif 350 < CO2_max < 400 and 120 < CO2_min < 150:
        points += 1
    else:
        points += 0
        
    if T_s_max < 275 and T_s_min> 257:
        points += 5
        print('full score for Ts !')

    elif 275 < T_s_max < 277 and 255 < T_s_min < 257:
        points += 2 
    elif 277 < T_s_max < 279 and 253 < T_s_min < 255:
        points += 1
    else:
        points += 0
    
    if m_max < 130 and m_min> -10:
        points += 5
        # print('full score for ice mass !')

    elif 130 < m_max < 160 and -30 < m_min < -10:
        points += 2 
        print('suboptimal score for ice mass!')
    elif 160 < m_max < 180 and -50 < m_min < -30:
        points += 1
        print('poor score for ice mass!')
    else:
        print('bad score for ice mass!')
        points += 0
    
    score = round((points/15 *10), 1)
    
    return score # between 0 and 10
    
def plot_func(time, CO2array, Tarray, marray, run_number):
        pl.plot(time, marray*Gt_to_SLE_conv, color='cyan', label = 'ice mass')
        pl.title('Run: '+ str(run_number), fontsize = 20)
        pl.xlabel('time [ka]', fontsize=14)
        pl.ylabel('ice mass [m SLE]', fontsize=14)
        pl.legend()
        pl.grid(True)
        pl.show()

        # # surface & upper ocean temp
        pl.plot(time,  Tarray, color='red', label = 'Surface temperature')
        pl.title('Run: ' + str(run_number), fontsize = 20)
        pl.xlabel('time [ka]', fontsize=14)
        pl.ylabel('temperature [K]', fontsize=14)
        pl.grid(True)
        pl.legend()
        pl.show()
        
        # CO_2 concentration
        pl.plot(time, CO2array, color='black')
        pl.title('Run: '+ str(run_number), fontsize = 20)
        pl.xlabel('time [ka]', fontsize=14)
        pl.ylabel(r'$CO_2$ concentration [ppm]', fontsize=14)
        pl.grid(True)
        pl.show()

        
        return


### main

def main_sim(coeffs_a, coeffs_b, coeffs_c, coeffs_d, coeffs_e, runtime, timestep):
    
    # initial values 
    T_o = 278. # K 
    T_s = 280. # K 
    m = 10**7 # Pg ice #TODO: find initial value based on paleoclimate data
    D = 1. # m 
    CO2 = 330.  # ppm
    D = 100 # m 
    I = Insol(0)
    CO2 = 330  # ppm
    
    # declare output arrays
    T_o_arr = np.zeros(int(time/dt+1))
    T_s_arr = np.zeros(int(time/dt +1))
    m_arr = np.zeros(int(time/dt +1))
    
    D_arr = np.zeros(int(time/dt +1))
    I_arr = np.zeros(int(time/dt +1))
    CO2_arr = np.zeros(int(time/dt +1))
    
    m_abl_arr = np.zeros(int(time/dt +1))
    acc_arr = np.zeros( int(time/dt +1))
    s_abl_arr = np.zeros(int(time/dt +1))
    
    comp1 = np.zeros(int(time/dt +1))
    comp2 = np.zeros(int(time/dt +1))
    comp3 = np.zeros(int(time/dt +1))

    
    for t in range(0, time+1, dt): 
        
        # first calculate new values for all model variables using the old values
        T_o_new = T_o + dTodt(T_s, T_o, coeffs_c)*dt
        I_new = Insol(t)
        D_new = D + dDdt(m, D, coeffs_d)*dt #TODO: adjust time scale?      
    
        CO2_new = f_CO2(CO2, T_o, coeffs_e )
        if CO2_new < CO2_min:
            CO2_new = CO2_min
        if CO2_new > CO2_max:
            CO2_new = CO2_max 
        
        dm_dt = dmdt(m, T_s, T_o, D, coeffs_a)
        
        m_new = m + dm_dt[0]*dt
        if m_new < 0:
            m_new = 0
    
        T_surface = T_surf(m, I, CO2, coeffs_b )
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
        # T_o_arr[t//10] = T_o
        T_s_arr[t//10] = T_s
        # D_arr[t//10] = D
        
        # I_arr[t//10] = I
        CO2_arr[t//10] = CO2
        m_arr[t//10] = m
        
        # m_abl_arr[t//10] = dm_dt[3]
        # s_abl_arr[t//10] = dm_dt[2]
        # acc_arr[ t//10] = dm_dt[1]
        
        # comp1[t//10] = T_surface[1]
        # comp2[t//10] = T_surface[2]
        # comp3[t//10] = T_surface[3]
        
    return  CO2_arr, T_s_arr, m_arr, coeffs_a, coeffs_b, coeffs_c, coeffs_d, coeffs_e



#%%

# define interval where we want to look for parameters 

# a_limits = [[a_0 / 10, a_0 *10], [a_1/10, a_1*10], [a_2/10, a_2*10], [a_3/10, a_3*10], [a_4/10, a_4*10], [a_5/10, a_5*10], [a_6/10, a_6 *10]]
a_limits = [[a_0 * 0.9, a_0 *1.1], [a_1 *0.9, a_1*1.1], [a_2*0.9, a_2*1.1], [a_3*0.9, a_3*1.1], [a_4*0.9, a_4*1.1], [a_5*0.9, a_5*1.1], [a_6*0.9, a_6 *1.1]]
b_limits = [[b_0 /100, b_0 * 100], [b_1/10, b_1*10], [b_2/10, b_2*10], [b_3/10, b_3*10]]
c_limits = [[c_0/10, c_0*10]]
d_limits = [[d_0/10, d_0*10],[d_1/10, d_1*10]]
e_limits = [[e_0/10, e_0*10], [e_1/10, e_1*10], [e_2/10, e_2*10]]

limits = []

limits.extend(a_limits)
limits.extend(b_limits)
limits.extend(c_limits)
limits.extend(d_limits)
limits.extend(e_limits)

limits = np.asarray(limits)
# use latin "hyper" cube sampling ( just a latin square for now)

sampling = LHS(xlimits=limits)
coefficients = sampling(n_runs)



results =  [[] for x in range(n_runs)]

# produce time axis
t_axis_f = np.arange(0,time+1,dt) #year, counting up from starting point
t_axis_r = np.zeros(len(t_axis_f)) #years ago
for j in range(len(t_axis_r)):
    t_axis_r[j] = t_axis_f[-(j+1)]

for run in range(n_runs): 
    
    run_simulation = main_sim(coefficients[run][:7], coefficients[run][7:11], coefficients[run][11], coefficients[run][12:14], coefficients[run][14:17], time, dt)
    # run_simulation = main_sim(a_coeffs, b_coeffs, c_0, d_coeffs, e_coeffs, time, dt)
    
    score = scorefunction(np.asarray(run_simulation[:3]))
    
    if score > 0:
        results[run].append( run_simulation )
        results[run].append(score)
        print('score: ', score, 'for run number: ' + str(run))
        plot_func(t_axis_f, results[run][0][0], results[run][0][1], results[run][0][2], run)

    else:
        results[run].append(np.NaN)
        results[run].append(np.NaN)
        




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
import pandas as pd
import csv


# simulation parameters
dt = 10 # yrs
time = int(1000 * 1e3) # yrs

n_runs = 2500 # amount of simulation runs out of which the top 5 is selected
n_iterations = 1000 # how many iterations of n_runs runs

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
CO2_min = 5
# CO2_min = e_2 
CO2_max = 1600.
# CO2_max = 600. 
e_coeffs = np.array((e_0, e_1,e_2))

#other
P_max = 3*10**-3 #Pg per km^2 per yr

# P_max = P_max * Gt_to_SLE_conv 

T_ref = 273.15
T_min = 228.15

# function that calculates dm/dt, input list of a coefficients as coef
def dmdt(m, T_s, T_o, D, coef):
    m_max = 100000 # MSLE # this value is here just to keep the model from running into overflows and crashing the loop 
    
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
                                               
    if m * Gt_to_SLE_conv > m_max:
        mass_change= 0
    else: 
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

def penaltyfunction(arrays, arrays2): 
    # assuming arrays is a tuple object with arrays for CO_2, T_s and  ice mass
    
    # exlude spinup time of 50 ka: 50*1e3 /10 = 5000 timesteps
    spinup = 5000
    
    CO2_max = np.max(arrays[0, spinup:])
    CO2_min = np.min(arrays[0, spinup:])
    
    T_s_max = np.max(arrays[1, spinup:])
    T_s_min = np.min(arrays[1, spinup:])
    
    m_max = np.max(arrays[2, spinup:]*Gt_to_SLE_conv)
    m_min = np.min(arrays[2, spinup:]*Gt_to_SLE_conv)
    
    # next for the later added arrays
    mabl_max = np.max(arrays2[0, spinup:]*Gt_to_SLE_conv)
    mabl_min = np.min(arrays2[0, spinup:]*Gt_to_SLE_conv)
    
    sabl_max = np.max(arrays2[1, spinup:]*Gt_to_SLE_conv)
    sabl_min = np.min(arrays2[1, spinup:]*Gt_to_SLE_conv)
    
    acc_max = np.max(arrays2[2, spinup:]*Gt_to_SLE_conv)
    acc_min = np.min(arrays2[2, spinup:]*Gt_to_SLE_conv)
    
    
    penalty = 0
    
    #CO2
    penalty += 0.01 * (CO2_max - 300)**2
    penalty += 0.01 * (CO2_min - 180)**2

    # T_s
    penalty += 0.01 * (T_s_max - 275)**2
    penalty += 0.01 * (T_s_min - 257)**2

    # m
    penalty += 0.01 * (m_max - 130)**2
    penalty += 0.01 * (m_min - 10)**2
    
    #added these 3 later
    # m abl 
    penalty += 0.01 * (mabl_max + 1*0.025)**2
    penalty += 0.01 * (mabl_min + 1*0.025)**2
    
    # s abl
    penalty += 0.01 * (mabl_max + 1*0.025)**2
    penalty += 0.01 * (mabl_min + 1*0.025)**2
    
    # acc
    penalty += 0.01 * (mabl_max - 3*0.025)**2
    penalty += 0.01 * (mabl_min - 3*0.025)**2
    
    # from a back of the envelope calculation: the mass balance must be about
    # 9025 Gt/10yrs = 0.025 m SLE/10yrs to accumulate realistic amounts of mass. Here we assume
    # a ratio of accumulation to ablation of 3:2 in total. This is a 
    # guess. 
    
    result = round((penalty), 1)
    
    return result # between 0 and 10
    
def plot_func(time, CO2array, Tarray, marray, run_number, score):
        pl.plot(time, marray*Gt_to_SLE_conv, color='cyan', label = 'ice mass')
        pl.title('Run: '+ str(run_number)+ ', score ' +str(score), fontsize = 20)
        pl.xlabel('time [ka]', fontsize=14)
        pl.ylabel('ice mass [m SLE]', fontsize=14)
        pl.legend()
        pl.grid(True)
        pl.show()

        # # surface & upper ocean temp
        pl.plot(time,  Tarray, color='red', label = 'Surface temperature')
        pl.title('Run: ' + str(run_number)+ ', score ' +str(score), fontsize = 20)
        pl.xlabel('time [ka]', fontsize=14)
        pl.ylabel('temperature [K]', fontsize=14)
        pl.grid(True)
        pl.legend()
        pl.show()
        
        # CO_2 concentration
        pl.plot(time, CO2array, color='black')
        pl.title('Run: '+ str(run_number)+ ', score ' +str(score), fontsize = 20)
        pl.xlabel('time [ka]', fontsize=14)
        pl.ylabel(r'$CO_2$ concentration [ppm]', fontsize=14)
        pl.grid(True)
        pl.show()

        
        return

def plot_func_mass_components(time, m, mabl, sabl, acc, run_number, score):
        
        pl.plot(time, m*Gt_to_SLE_conv, color='green', label = 'ice mass')

        # # marine ablation 
        pl.plot(time/1000,  mabl * Gt_to_SLE_conv, color='orange', label = 'marine ablation')
        # surface ablation
        pl.plot(time/1000, sabl * Gt_to_SLE_conv , color='red', label = 'surface ablation')
        # accumulation
        pl.plot(time/1000, acc * Gt_to_SLE_conv, color='blue', label = 'accumulation')

        pl.legend()

        pl.title('Run: '+ str(run_number)+ ', score' +str(score), fontsize = 20)
        pl.xlabel('time [ka]', fontsize=14)
        pl.ylabel('ice mass change', fontsize=14)
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
        
        m_abl_arr[t//10] = dm_dt[3]
        s_abl_arr[t//10] = dm_dt[2]
        acc_arr[ t//10] = dm_dt[1]
        
        # comp1[t//10] = T_surface[1]
        # comp2[t//10] = T_surface[2]
        # comp3[t//10] = T_surface[3]
        
    return  CO2_arr, T_s_arr, m_arr, coeffs_a, coeffs_b, coeffs_c, coeffs_d, coeffs_e, m_abl_arr, s_abl_arr, acc_arr

#%%

# define interval where we want to look for parameters, right now +- 10% 

### a_limits = [[a_0 / 10, a_0 *10], [a_1/10, a_1*10], [a_2/10, a_2*10], [a_3/10, a_3*10], [a_4/10, a_4*10], [a_5/10, a_5*10], [a_6/10, a_6 *10]]
# a_limits = [[a_0 * 0.7, a_0 * 1.3], [a_1 * 0.7, a_1 * 1.3], [a_2 * 0.7, a_2 * 1.3], [a_3 * 0.7, a_3 * 1.3], [a_4 * 0.7, a_4 * 1.3], [a_5 * 0.7, a_5 * 1.3], [a_6 * 0.7, a_6 * 1.3]]
# b_limits = [[b_0 * 0.7, b_0 * 1.3], [b_1 * 0.7, b_1 * 1.3], [b_2 * 0.7, b_2 * 1.3], [b_3 * 0.7, b_3 * 1.3]]
# c_limits = [[c_0 * 0.7, c_0 * 1.3]]
# d_limits = [[d_0 * 0.7, d_0 * 1.3], [d_1 * 0.7, d_1 * 1.3]]
# e_limits = [[e_0 * 0.7, e_0 * 1.3], [e_1 * 0.7, e_1 * 1.3], [e_2 * 0.7, e_2 * 1.3]]

# +- 10 % 
a_limits = [[a_0 * 0.9, a_0 * 1.1], [a_1 * 0.9, a_1 * 1.1], [a_2 * 0.9, a_2 * 1.1], [a_3 * 0.9, a_3 * 1.1], [a_4 * 0.9, a_4 * 1.1], [a_5 * 0.9, a_5 * 1.1], [a_6 * 0.9, a_6 * 1.1]]
b_limits = [[b_0 * 0.9, b_0 * 1.1], [b_1 * 0.9, b_1 * 1.1], [b_2 * 0.9, b_2 * 1.1], [b_3 * 0.9, b_3 * 1.1]]
c_limits = [[c_0 * 0.9, c_0 * 1.1]]
d_limits = [[d_0 * 0.9, d_0 * 1.1], [d_1 * 0.9, d_1 * 1.1]]
e_limits = [[e_0 * 0.9, e_0 * 1.1], [e_1 * 0.9, e_1 * 1.1], [e_2 * 0.9, e_2 * 1.1]]


    # produce time axis
t_axis_f = np.arange(0,time+1,dt) #year, counting up from starting point
t_axis_r = np.zeros(len(t_axis_f)) #years ago
for j in range(len(t_axis_r)):
    t_axis_r[j] = t_axis_f[-(j+1)]

def make_limits_from_minmax(array):
    array = np.asarray(array)
    # rows = array.shape[0]
    columns = array.shape[1]
    
    lims = [[] for x in range(columns)]
    for col in range(columns):
        lims[col].append(np.min(array[:,col]))
        lims[col].append(np.max(array[:,col]))       
    
    return lims

def parameter_search(a_lims, b_lims, c_lims, d_lims, e_lims, runN, time, dt):
    limits = []
    
    limits.extend(a_lims)
    limits.extend(b_lims)
    limits.extend(c_lims)
    limits.extend(d_lims)
    limits.extend(e_lims)
    
    limits = np.asarray(limits)
    
    ## force minimum 1% search range
    for i in range(len(limits)):
        low_bound = limits[i][0]
        high_bound = limits[i][1]
        
        diff = np.abs(high_bound - low_bound)
                
        if diff < 0.01 * low_bound or diff < 0.01*high_bound:
            new_low = 0.999 * low_bound
            new_high = 1.001 * high_bound
            
            limits[i][0] = new_low
            limits[i][1] = new_high
            
        else: 
            pass
###    use latin hyper cube sampling
    
    sampling = LHS(xlimits=limits)
    coefficients = sampling(runN)
    
    results =  [[] for x in range(runN)]
    
    run_and_score = np.zeros((2,runN))
    for run in range(runN): 
        
        run_simulation = main_sim(coefficients[run][:7], coefficients[run][7:11], coefficients[run][11], coefficients[run][12:14], coefficients[run][14:17], time, dt)
        # run_simulation = main_sim(a_coeffs, b_coeffs, c_0, d_coeffs, e_coeffs, time, dt)
        
        penalty = penaltyfunction(np.asarray(run_simulation[:3]), np.asarray(run_simulation[8:12]))
            
        results[run].extend( run_simulation)
        results[run].append(penalty)
        if run % 25 == 0:
            print(f'penalty points: {penalty} for run number {run}')
        # print(f'penalty points: {penalty} for run number {run}')
        # plot_func(t_axis_f, results[run][0], results[run][1], results[run][2], run)
    
        run_and_score[0,run] = int(run)
        run_and_score[1,run]= penalty 
        
    if penalty == 0: 
            print(f'*!!!* Run number {run} has completed without incurring any penalties*!!!*')
                  
    df_results = pd.DataFrame(np.transpose(run_and_score), columns=['run', 'score'])
    df_results = df_results.sort_values('score')
    
    top_5_runs = df_results.iloc[:5,0]
    top_5_data = []
    topscore = df_results.iloc[0,1]
    
    best_a = np.zeros((5, len(a_coeffs)))
    best_b = np.zeros((5, len(b_coeffs)))
    best_c = np.zeros((5, 1))
    best_d = np.zeros((5, len(d_coeffs)))
    best_e = np.zeros((5, len(e_coeffs)))
    
    for i in enumerate(top_5_runs):    
        top_5_data.append(results[int(i[1])]) # not sure if we need this data seperately
        
        best_a[i[0],:] = results[int(i[1])][3] 
        best_b[i[0],:] = results[int(i[1])][4]
        best_c[i[0],:] = results[int(i[1])][5]
        best_d[i[0],:] = results[int(i[1])][6]
        best_e[i[0],:] = results[int(i[1])][7]

    new_a_range = make_limits_from_minmax(best_a)
    new_b_range = make_limits_from_minmax(best_b)
    new_c_range = make_limits_from_minmax(best_c)
    new_d_range = make_limits_from_minmax(best_d)
    new_e_range = make_limits_from_minmax(best_e)

    return new_a_range, new_b_range, new_c_range, new_d_range, new_e_range, top_5_data, topscore

# calculate starting point from best guess made by hand
firstguess = parameter_search(a_limits, b_limits, c_limits, d_limits, e_limits, 1, time, dt)

a = firstguess[0]
b = firstguess[1]
c = firstguess[2]
d = firstguess[3]
e = firstguess[4]

# create file, set headers
with open('parametersearchresults2500_1000ka2.csv', mode='w') as csvfile:
    header_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    header_writer.writerow(f'nrun is {n_runs}, iterations: {n_iterations}')
    header_writer.writerow(['a0', 'a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'b0', 'b1', 'b2', 'b3', 'c0', 'd0', 'd1', 'e0', 'e1', 'e2', 'score'])
    
for iterations in range(n_iterations):
    print(f'ITERATION #: {iterations}')
    simulation = parameter_search(a, b, c, d, e, n_runs, time, dt )
    
    a = simulation[0]
    b = simulation[1]
    c = simulation[2]
    d = simulation[3]
    e = simulation[4]    
    
    output = []
    output.append(simulation[:5])
    with open('parametersearchresults2500_1000ka2.csv', mode='a') as csvfile:
        limits_writer = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        output.append(simulation[6])
        limits_writer.writerow(output)
    print(f'This iterations topscore is: {simulation[6]}')

#%%


convergence =[[-5.195642692826793e-08, -5.0591248467859786e-08], [6.975611737995667e-06, 6.988602729324854e-06], [3.268700120170882, 3.271495626848382], [34730.27770011655, 34945.454950208055], [0.00020349036753745237, 0.00020603577622953424], [3.607396443930929e-07, 3.8280524307642596e-07], [7.160152933980525e-08, 7.167329255126716e-08]], [[89.44808101316544, 89.6534142840577], [0.010387217564774783, 0.010669446552885046], [4.104564417647632e-10, 4.129099836252062e-10], [5.137865736133415, 5.203229220622828]], [[1.589468846108354e-07, 1.5941961158648341e-07]], [[2.5473936146013737e-09, 2.5790274803604174e-09], [7879.951286582043, 7915.4188516570675]], [[97.5698243255616, 100.67667004171584], [8.996539735677258, 8.999473377760545], [90.77755592526381, 90.81367328584854]]
conv_score = 358.7


aa = np.mean(np.asarray(convergence[0]), axis=1)
bb = np.mean(np.asarray(convergence[1]), axis=1)
cc = np.mean(np.asarray(convergence[2]), axis=1)
dd = np.mean(np.asarray(convergence[3]), axis=1)
ee = np.mean(np.asarray(convergence[4]), axis=1)
              
run_conv = main_sim(aa, bb, cc, dd, ee, time, dt)

plot_func(t_axis_f, run_conv[0], run_conv[1], run_conv[2], 0, conv_score )
# plot_func_mass_components(t_axis_f, run_conv[2], run_conv[8], run_conv[9], run_conv[10], 0)

### loop 2 
#### iterate searching 


# n_2_runs = 1   
# results2  =  [[] for x in range(n_runs)]

# run_and_score2 = np.zeros((2,n_runs))
# for run in range(n_2_runs): 
#     run_simulation2 = main_sim(new_a, new_b, new_c, new_d, new_e, time, dt)
#     penalty2 = penaltyfunction(np.asarray(run_simulation2[:3]))

#     results2[run].extend( run_simulation2)
#     results2[run].append(penalty2)
    
#     print(f'penalty points: {penalty2} for run number {run}')
#     plot_func(t_axis_f, results2[run][0], results2[run][1], results2[run][2], run)


#     run_and_score2[0,run] = int(run)
#     run_and_score2[1,run]= penalty2     

# #%%
    
# m_abl = run_simulation2[8]
# s_abl = run_simulation2[9]
# acc =  run_simulation2[10]

m_arr= run_conv[2]
m_abl = run_conv[8]
s_abl = run_conv[9]
acc = run_conv[10]
    
# pl.plot(t_axis_f, m_arr*Gt_to_SLE_conv, color='cyan', label = 'ice mass')
pl.plot(t_axis_f/1000,m_abl*Gt_to_SLE_conv, color='red', label = 'm abl ')
pl.plot(t_axis_f/1000, s_abl*Gt_to_SLE_conv, color='yellow', label = 's abl')
pl.plot(t_axis_f/1000,acc*Gt_to_SLE_conv, color='green', label = 'acc')
pl.title('? run, score: '+str(conv_score), fontsize = 20)
pl.xlabel('time [ka]', fontsize=14)
pl.ylabel('ice mass change [m SLE / 10 yr]', fontsize=14)
pl.legend()
pl.grid(True)
pl.show()




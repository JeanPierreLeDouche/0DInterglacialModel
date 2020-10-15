# -*- coding: utf-8 -*-
"""
Created on Thu Oct  8 11:47:29 2020

@author: Gebruiker
"""

import numpy as np

### Impressive arithmatic

a = 2
b = 3 
c = a * b 

### cutting edge modelling: 

def tempfunc(insolation, CO2_conc, mass_ice, c_1, c_2, c_3):
    T = c_1 * insolation + c_2 * CO2_conc + c_3 * mass_ice
    return T


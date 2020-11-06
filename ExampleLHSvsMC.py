# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 11:15:05 2020

@author: Gebruiker
"""

import numpy as np
import matplotlib.pyplot as plt

from smt.sampling_methods import LHS
import random as r

# use latin "hyper" cube sampling ( just a latin square for now)
xlimits = np.array([[0.0, 4.0], [0.0, 3.0]])
sampling = LHS(xlimits=xlimits)

num = 50
x = sampling(num)

# also monte carlo sampling
L_x_mc= []
L_y_mc=[]
for i in range(50):
    x_mc = r.uniform(0, 4)
    y_mc = r.uniform(0, 3)
    L_x_mc.append(x_mc)
    L_y_mc.append(y_mc)
#%%

# plt.title('LHS sampling vs MC sampling example', fontsize = 25)
plt.plot(x[:, 0], x[:, 1], "o", color = 'b', label = 'LHS')
plt.plot(L_x_mc, L_y_mc, "x", color = 'r', label = 'MC')
plt.xlabel("x", fontsize = 20)
plt.xticks(fontsize=15)
plt.ylabel("y", fontsize = 20)
plt.yticks(fontsize = 15)
plt.legend(fontsize = 15)
plt.show()
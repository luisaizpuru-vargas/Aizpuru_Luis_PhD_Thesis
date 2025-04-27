#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 24 12:48:46 2025

@author: luis
"""

#%%


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

x = [0,1,2,3]


# r = 0.3
y_EOC = [3, 64, 4305, 15507 ]

y_BBS = [7, 941, 12466, 17484]

y_HS = [4, 92, 3796 , 15700 ]

# Plot the curves
plt.plot(x, y_EOC, label = 'Early Onset Cardiomyopathy' )    # Curve 1
plt.plot(x, y_BBS, label = 'Beardet Biedl Syndrome' ) 
plt.plot(x, y_HS, label = 'Hirschsprung Disease' ) 


# Add data points
plt.scatter(x, y_EOC, s=50, marker='o')
plt.scatter(x, y_BBS, s=50, marker='o')
plt.scatter(x, y_HS, s=50, marker='o')


# Set y-axis to logarithmic scale
plt.yscale('log')
plt.ylim(-100, 70000)

# Customize x-axis ticks
plt.xticks(ticks=[0, 1, 2, 3], labels=["0th order", "1st order", "2nd order", "3rd order"])

# Customize the plot
plt.title('Unique Peripherals for 3 Oligogenic Diseases at |r|>0.3', fontsize = 12)
plt.xlabel('Order of Peripheral', fontsize = 12)
plt.ylabel('Cumulative Number of Unique Peripherals', fontsize = 10.5)
plt.legend()
plt.grid(True, which="major", ls="--")  # Grid for both major and minor ticks

# Display the plot
plt.show()
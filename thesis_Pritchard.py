#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 28 10:37:56 2025

@author: luis
"""
#### Thesis 1 ############


import numpy as np
import matplotlib.pyplot as plt

# Set parameters
a = 1000
rate = 1/a        # Rate parameter for the exponential distribution
num_samples = 10000 # Number of samples to generate
num_bins = 1000   # Number of bins in the histogram

# Generate random samples from an exponential distribution
samples = np.random.exponential(scale=5/rate, size=num_samples)

# Create the histogram
plt.figure(figsize=(10, 6))  # Set figure size for better visibility
plt.hist(samples, bins=num_bins, density=True, alpha=0.7, color='blue')

# Remove ticks on x and y axes
plt.xticks([2500], ['a'], fontsize = 18)
plt.yticks([])

# Add labels and title (unitless)
plt.xlabel('Polymorphic Loci', fontsize=18)  # No units specified
plt.ylabel('H^2', fontsize=18)              # No units specified
plt.title('Histogram of H^2 as a Function of Polymorphic Loci', fontsize=18)

#plt.yscale("log")
plt.xlim(0, 30000)

# Add a vertical line at x=2500
plt.axvline(x=2500, color='red', linestyle='--', label='Threshold for generally considered Loci', linewidth = 2.5)

# Add the legend
plt.legend(fontsize=16)

# Display the plot
plt.show()


#%%


##### Thesis 2 ###########3

import numpy as np
import matplotlib.pyplot as plt

num_bins = 23
points_per_bin = 10000

# Generate random heights for each bin between 0 and 0.07
h = np.random.uniform(0.038, 0.046, num_bins)

data = []
for i in range(1, num_bins + 1):
    low = i - 0.5
    high = i + 0.5
    bin_data = np.random.uniform(low, high, points_per_bin)
    data.extend(bin_data)
data = np.array(data)

# Create weights: each point's weight is h[i] / points_per_bin for its bin
weights = np.repeat(h / points_per_bin, points_per_bin)

# Define bin edges from 0.5 to 23.5
bin_edges = np.arange(0.5, 24.5, 1)

plt.figure(figsize=(10, 6))
plt.hist(data, bins=bin_edges, weights=weights,  color='blue', alpha=0.7, rwidth = 0.7)
plt.xlabel('Chromosome', fontsize=18)
plt.ylabel('Heritability', fontsize =18)
plt.title('Histogram with Equal Contribution for H^2 per chromosome', fontsize =18)
plt.ylim(0,0.07)

# Set x-ticks at 1 to 23
plt.xticks(np.arange(1, 24, 1), [str(i) for i in range(1, 24)])

# Add a horizontal line at x=0.04
plt.axhline(y=0.04, color='red', linestyle='--', label='Average H^2 contribution per Chromosome', linewidth = 2.5)
# Add the legend
plt.legend(fontsize=16)


plt.yticks([0.045], ['(H^2)/23'], fontsize = 14)
plt.setp(plt.gca().get_yticklabels(), rotation=90)


plt.show()




#%%

## Thesis 3 ####

import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(0)

# Generate sample data for two groups
data1 = np.random.normal(0.58, 0.08, 100)  # Mean=0, Std=1, 100 samples
data2 = np.random.normal(0.42, 0.08, 100)  # Mean=1, Std=1.5, 100 samples

# Define positions for the box plots on the x-axis
positions = [1, 2]

# Create the box plots
plt.figure(figsize=(6, 4))  # Set figure size (optional)
plt.boxplot([data1, data2], positions=positions)


# Set major ticks every 0.1 from 0 to 1 on the y-axis
major_ticks = np.arange(0, 1.1, 0.1)
plt.yticks(major_ticks)

# Hide all y-tick labels
plt.gca().tick_params(axis='y', labelleft=False)


# Set x-ticks and labels
plt.xticks(positions, ['Relevant Genes', 'Random Genes'], fontsize = 13)


# Add labels and title
plt.xlabel('Type of Genes', fontsize=13)
plt.ylabel('H^2 per SNP', fontsize =13)
plt.title('Box Plots for H^2 contribution for relevant vs random genes', fontsize =12)
plt.ylim(0,1)

# Add a grid for better readability (optional)
plt.grid(True, axis='y')

# Display the plot
plt.show()




#%%


## Thesis 4 ####

import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(0)

# Generate sample data for two groups
data1 = np.random.normal(0.58, 0.06, 100) / 100  # Mean=0, Std=1, 100 samples
data2 = np.random.normal(0.42, 0.06, 100) / 100  # Mean=1, Std=1.5, 100 samples

# Define positions for the box plots on the x-axis
positions = [1, 2]

# Create the box plots
plt.figure(figsize=(6, 4))  # Set figure size (optional)
plt.boxplot([data1, data2], positions=positions)

# Set x-ticks and labels
plt.xticks(positions, ['Tissue Specific ', 'Broadly Active'], fontsize = 13)

# Add labels and title
plt.xlabel('Type of Regulatory Elements', fontsize=13)
plt.ylabel('H^2 per SNP', fontsize =13)
plt.title('Box Plots for H^2 per SNP in Regulatory Elements', fontsize =14)
plt.ylim(0.002,0.0075)


# Set major ticks every 0.1 from 0 to 1 on the y-axis
major_ticks = np.arange(0.002, 0.008, 0.0005)
plt.yticks(major_ticks)

# Hide all y-tick labels
plt.gca().tick_params(axis='y', labelleft=False)

# Add a grid for better readability (optional)
plt.grid(True, axis='y')

# Display the plot
plt.show()


#%%

### Thesis 5 #####


import matplotlib.pyplot as plt
import numpy as np

# Define the data
categories = ['Complex Phenotype', 'Mendellian Phenotype']  # Labels for the two ticks
group1_values = [.9, .1]  # Heights for Group 1 bars
group2_values = [.1, .9]  # Heights for Group 2 bars

# Set parameters
bar_width = 0.35  # Width of each bar
category_positions = np.arange(len(categories))  # Positions for the ticks (0, 1)

# Calculate bar positions
positions_group1 = category_positions - bar_width / 2  # Offset Group 1 bars to the left
positions_group2 = category_positions + bar_width / 2  # Offset Group 2 bars to the right

# Create the plot
plt.figure(figsize=(6, 4))  # Set figure size
plt.bar(positions_group1, group1_values, width=bar_width, label='Common Variants', color='blue')
plt.bar(positions_group2, group2_values, width=bar_width, label='Rare Variants', color='orange')

# Customize the plot
plt.xticks(category_positions, categories)  # Set ticks at the center of each group
plt.xlabel('Type of Phenotype', fontsize = 12)  # Label for x-axis
plt.ylabel('H^2', fontsize = 12)  # Label for y-axis
plt.title('Heritabilty by Type of Variant', fontsize = 12)  # Title
plt.legend()  # Add legend to distinguish groups


# Set major ticks every 0.1 from 0 to 1 on the y-axis
major_ticks = np.arange(0, 1, 0.1)
plt.yticks(major_ticks)

# Hide all y-tick labels
plt.gca().tick_params(axis='y', labelleft=False)

# Add a grid for better readability (optional)
plt.grid(True, axis='y')

# Display the plot
plt.show()


#%%
##### Thesis 6 ######



import numpy as np
import matplotlib.pyplot as plt

# Set random seed for reproducibility
np.random.seed(0)

# Generate sample data for two groups
data1 = np.random.normal(0.35, 0.05, 100)   # Mean=0, Std=1, 100 samples
data2 = np.random.normal(0.65, 0.05, 100)  # Mean=1, Std=1.5, 100 samples

# Define positions for the box plots on the x-axis
positions = [1, 2]

# Create the box plots
plt.figure(figsize=(6, 4))  # Set figure size (optional)
plt.boxplot([data1, data2], positions=positions)

# Set x-ticks and labels
plt.xticks(positions, ['Coding Variants', 'Non-coding Variants'], fontsize = 13)

# Add labels and title
plt.xlabel('Type of Variant', fontsize=13)
plt.ylabel('H^2 ', fontsize =13)
plt.title('Box Plots for H^2 for type of coding variants', fontsize =14)
#plt.ylim(0,1)


# Set major ticks every 0.1 from 0 to 1 on the y-axis
major_ticks = np.arange(0, 1, 0.1)
plt.yticks(major_ticks)

# Hide all y-tick labels
plt.gca().tick_params(axis='y', labelleft=False)

# Add a grid for better readability (optional)
plt.grid(True, axis='y')

# Display the plot
plt.show()












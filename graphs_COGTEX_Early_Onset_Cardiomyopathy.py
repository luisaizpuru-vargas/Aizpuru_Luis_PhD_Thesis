#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 12 14:25:05 2025

@author: luis
"""


import os
import ctypes
import subprocess
import numpy as np

    
#%%

#### Importar el .RData ####

## Importar Librerias ###

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri

## enables automatic conversion of R data frames to pandas DataFrames
pandas2ri.activate()

file_path = "/home/luis/Desktop/ITESM_PhD/Semesters/2ndo_Semestre/Defensa_Propuesta/Data/Z_score_pearson_min.RData"

# Load .RData File
r = robjects.r
r.load(file_path)

## List Loaded Objects
loaded_objects = list(r("ls()"))
print("Loaded objects:", loaded_objects)


# Import object m
m = robjects.r['m']
m_python = np.array(m)
print(len(m_python))

#%% Summary Statistics

# Check its class
obj_class = robjects.r('class(m)')[0]
print(f"The class of 'm' is: {obj_class}")

print("Length of m:", len(m_python))

print("Min:", np.min(m_python))
print("Max:", np.max(m_python))
print("Mean:", np.mean(m_python))
print("Median:", np.median(m_python))
print("Standard deviation:", np.std(m_python))


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


# Calculate statistics
min_val = np.min(m_python)/100
max_val = np.max(m_python)/100
mean_val = np.mean(m_python)/100
median_val = np.median(m_python)/100
std_val = np.std(m_python)/100

# Create the histogram
plt.figure(figsize=(10, 6))  # Set figure size for better visibility
plt.hist(m_python/100, bins=100, alpha=0.7, color='blue', edgecolor='black')

# Add vertical lines for mean and median
plt.axvline(mean_val, color='red', linestyle='dashed', linewidth=3)
plt.axvline(median_val, color='green', linestyle='dashed', linewidth=3)


# Create custom legend entries
legend_elements = [
    Line2D([0], [0], color='red', linestyle='dashed', linewidth=2, label=f'Mean: {mean_val:.2f}'),
    Line2D([0], [0], color='green', linestyle='dashed', linewidth=2, label=f'Median: {median_val:.2f}'),
    Line2D([0], [0], color='white', label=f'Min: {min_val:.2f}'),
    Line2D([0], [0], color='white', label=f'Max: {max_val:.2f}'),
    Line2D([0], [0], color='white', label=f'Std: {std_val:.2f}'),
    Line2D([0], [0], color='white', label=f'# of elements: {560000000:.2e}')
]

# Add the legend
plt.legend(handles=legend_elements, loc='upper right')

# Add titles and labels
plt.title('Histogram for correlation values in CoGTEx v1.0')
plt.xlabel('Value')
plt.ylabel('Frequency')


## Log y axis
plt.yscale('log')

# Display the plot
plt.show()
   


#%%

######################################################################################
######################################################################################
######################################################################################
####################### PART #2 ######################################################




###### 1st Order Peripheral ###################3



### Query for associations ####

import matplotlib.pyplot as plt
import pandas as pd

# Step 1: Read headers from "geneinfoHeaders.txt"
with open("geneinfoHeaders.txt", "r") as f:
    headers = [line.strip() for line in f]

# Step 2: Read "geneinfo.txt" into a DataFrame and set column names
x = pd.read_csv("geneinfo.txt", sep="\t", header=None, names=headers)

# Step 3: Get the number of rows in x
n = x.shape[0]
#print(n)


#Step 4
# Import object m
m = robjects.r['m']
m_python = np.array(m)

# Step 5: Define the genes of interest
query = ["MKL2", "MYH7", "NKX2-5"]

# Step 6: Find the 1-based indices of query genes in 'Gene Symbol' column

cogtexIdxs = [int(x[x['Gene Symbol'] == gene].index[0]) + 1 for gene in query]



def getIndexPairs_numpy(idx, m, n):
    pidx = idx - 1  # Convert 1-based idx to 0-based
    associations = np.zeros(n)
    for j in range(n):
        if j == pidx:
            associations[j] = 0.0  # Distance to self is 0
        else:
            a, b = min(pidx, j), max(pidx, j)
            pair_idx = a * (n - 1) - (a * (a - 1)) // 2 + (b - a - 1)
            associations[j] = m[pair_idx]
    return associations


# Step 7: Extract and scale associations
estimates = []
for idx in cogtexIdxs:
    associations = getIndexPairs_numpy(idx, m_python, n)
    scaled_associations = associations / 100  # Scale as in the original script
    estimates.append(scaled_associations)    
    

idx = cogtexIdxs[0]
associations = getIndexPairs_numpy(idx, m_python, n)
#print(f"Raw associations: {associations}")
#print(f"Scaled associations: {associations / 100}")
#print(m_python)


# Step 8: Convert to pandas Series with gene symbols
gene_symbols = x['Gene Symbol'].tolist()
estimates_py = [pd.Series(est, index=gene_symbols) for est in estimates]


# Step 9: Create a dictionary with query genes as keys
estimates_dict = dict(zip(query, estimates_py))


#%%%%

###### Iteration 3 #############################
#########33 Counting of 1st order peripheral genes vs threshold ####################3

import numpy as np
from collections import Counter
import matplotlib.pyplot as plt

# Define parameters
thresholds = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9]
use_absolute = True  # Set to True to use absolute correlation values

# Lists to store data for plotting
threshold_values = []
unique_counts = []
total_counts = []
difference_counts = []
percentage_differences = []

# Loop over each threshold value
for threshold in thresholds:
    total_count = 0
    unique_coexpressed = set()
    
    # Analyze each queried gene
    for gene, est in estimates_dict.items():
        # Define condition based on whether to use absolute values
        if use_absolute:
            condition = (est.abs() > threshold) & (est.index != gene)
        else:
            condition = (est > threshold) & (est.index != gene)
        
        # Count co-expressed genes for this gene
        count = condition.sum()
        total_count += count
        
        # Update the set of unique co-expressed genes
        coexpressed_genes = est[condition].index
        unique_coexpressed.update(coexpressed_genes)
    
    # Calculate unique_count for this threshold
    unique_count = len(unique_coexpressed)
    
    # Store the threshold and the count of unique co-expressed genes
    threshold_values.append(threshold)
    unique_counts.append(unique_count)
    total_counts.append(total_count)
    difference_counts.append(total_count - unique_count)
    
    # Calculate percentage difference using the current unique_count
    if unique_count > 0:
        percentage_diff = ((total_count - unique_count) / unique_count) * 100
    else:
        percentage_diff = 0  # Avoid division by zero
    percentage_differences.append(percentage_diff)
    
    # Output summary results for this threshold
    print(f"\nThreshold: {threshold}")
    print(f"Total unique co-expressed genes with {'|r|' if use_absolute else 'r'} > {threshold} across all queried genes: {unique_count}")
    print(f"Total co-expressed genes (with duplicates) with {'|r|' if use_absolute else 'r'} > {threshold}: {total_count}")


# Create subplots: 2 rows, 1 column
fig, ax = plt.subplots(2, 1, figsize=(10, 10))  # figsize controls the width and height

# Plot 1: Unique Counts (top plot)
ax[0].plot(thresholds, unique_counts, marker='o', linestyle='-', color='blue', label = 'Unique Counts')
ax[0].set_title('1st Order Unique Counts vs. Threshold for Early-Onset Cardiomyopathy', fontsize = 16)
ax[0].set_xlabel('Threshold', fontsize = 14)
ax[0].set_ylabel('Unique Counts', fontsize = 14)
ax[0].grid(True)  # Add a grid for readability
ax[0].legend(loc='upper right', fontsize = 14)



# Plot 2: Percentage Differences (bottom plot)
ax[1].plot(thresholds, percentage_differences, marker='o', linestyle='-', color='green', label = "Percentage Difference")
ax[1].set_title('Percentage Difference (Total vs Unique)  vs. Threshold for 1st order Unique Counts' , fontsize = 16)
ax[1].set_xlabel('Threshold', fontsize = 14)
ax[1].set_ylabel('Percentage Difference (%)', fontsize = 14)
ax[1].grid(True)  # Add a grid for readability
ax[1].legend(loc='upper right', fontsize = 14)

# Adjust layout to prevent overlap and display the plots
plt.tight_layout()
plt.show()

#%%%


############################## 1st, 2nd and 3rd Order Peripherals #######################33

# Parameters
threshold = 0.9  # Correlation threshold (adjust as needed)
queried_genes =  set(["MKL2", "MYH7", "NKX2-5"])  # Example queried genes query = ["MKL2", "MYH7", "NKX2-5"]
use_absolute = True  # Set to True to use absolute correlation values


print(len(queried_genes))

#################################################
# Analyze each queried gene
### 1st order peripheral #####
########################################3
##########################################3
##############################################
 
first_order_genes = set()
for gene, est in estimates_dict.items():
     # Define condition based on whether to use absolute values
     if use_absolute:
         condition = (est.abs() > threshold) & (est.index != gene)
     else:
         condition = (est > threshold) & (est.index != gene)
     
     # Count co-expressed genes for this gene
     count = condition.sum()
     total_count += count
     
     # Update the set of unique co-expressed genes
     coexpressed_genes = est[condition].index
     first_order_candidates = set(coexpressed_genes) - queried_genes
     
     first_order_genes.update(first_order_candidates)

unique_count_1st_order = len(first_order_genes)
 
# Store the threshold and the count of unique co-expressed genes


# Output summary results for this threshold
print(f"\nThreshold: {threshold}")
print(f"Total unique 1st order co-expressed genes with {'|r|' if use_absolute else 'r'} > {threshold} across all queried genes: {unique_count_1st_order}")

#print(first_order_genes)


#%%
import pandas as pd
import numpy as np
from tqdm import tqdm

#%%
####################################################3
##################### 2nd order ##################3
####################################################

# Step 1: Read headers from "geneinfoHeaders.txt"
with open("geneinfoHeaders.txt", "r") as f:
    headers = [line.strip() for line in f]

# Step 2: Read "geneinfo.txt" into a DataFrame and set column names
x = pd.read_csv("geneinfo.txt", sep="\t", header=None, names=headers)

# Step 3: Get the number of rows in x
n = x.shape[0]
#print(n)


#Step 4
# Import object m
m = robjects.r['m']
m_python = np.array(m)

# Step 5: Define the genes of interest
query_1storder = first_order_genes
query_list_1storder = list(query_1storder)
#print(type(query_list_1storder))

# Step 6: Find the 1-based indices of query genes in 'Gene Symbol' column

cogtexIdxs = [int(x[x['Gene Symbol'] == gene].index[0]) + 1 for gene in query_list_1storder]

def getIndexPairs_numpy(idx, m, n):
    pidx = idx - 1  # Convert 1-based idx to 0-based
    associations = np.zeros(n)
    for j in range(n):
        if j == pidx:
            associations[j] = 0.0  # Distance to self is 0
        else:
            a, b = min(pidx, j), max(pidx, j)
            pair_idx = a * (n - 1) - (a * (a - 1)) // 2 + (b - a - 1)
            associations[j] = m[pair_idx]
    return associations


# Step 7: Extract and scale associations
estimates = []
for idx in tqdm(cogtexIdxs, desc="Extracting Associations"):
    associations = getIndexPairs_numpy(idx, m_python, n)
    scaled_associations = associations / 100  # Scale as in the original script
    estimates.append(scaled_associations)
    
    
idx = cogtexIdxs[0]

#print(idx)


#associations = getIndexPairs_numpy(idx, m_python, n)
#print(f"Raw associations: {associations}")
#print(f"Scaled associations: {associations / 100}")
#print(m_python)    

# Step 8: Convert to pandas Series with gene symbols
gene_symbols = x['Gene Symbol'].tolist()
estimates_py = [pd.Series(est, index=gene_symbols) for est in estimates]


# Step 9: Create a dictionary with query genes as keys
estimates_dict_1st = dict(zip(query_list_1storder, estimates_py))
    

queried_genes =  set(["MKL2", "MYH7", "NKX2-5"])
second_order_genes = set()


for gene in tqdm(first_order_genes, desc="Processing 2nd-Order Genes", leave=False):
    if gene in estimates_dict_1st:
        est = estimates_dict_1st[gene]
        if use_absolute:
            condition = (est.abs() > threshold) & (est.index != gene)
        else:
            condition = (est > threshold) & (est.index != gene)
        
        coexpressed_genes = est[condition].index
        second_order_candidates = set(coexpressed_genes) - queried_genes - first_order_genes
        second_order_genes.update(second_order_candidates)

unique_count_2nd_order = len(second_order_genes)


    
    
print(f"Number of unique 2nd-order peripheral genes: {unique_count_2nd_order}")



#%%%

####################################################3
##################### 3rd order ##################3
####################################################




# Step 1: Read headers from "geneinfoHeaders.txt"
with open("geneinfoHeaders.txt", "r") as f:
    headers = [line.strip() for line in f]

# Step 2: Read "geneinfo.txt" into a DataFrame and set column names
x = pd.read_csv("geneinfo.txt", sep="\t", header=None, names=headers)

# Step 3: Get the number of rows in x
n = x.shape[0]
#print(n)


#Step 4
# Import object m
m = robjects.r['m']
m_python = np.array(m)

# Step 5: Define the genes of interest
query_2ndorder = second_order_genes
query_list_2ndorder = list(query_2ndorder)
#print(type(query_list_2ndorder))

# Step 6: Find the 1-based indices of query genes in 'Gene Symbol' column

cogtexIdxs = [int(x[x['Gene Symbol'] == gene].index[0]) + 1 for gene in query_list_2ndorder]

def getIndexPairs_numpy(idx, m, n):
    pidx = idx - 1  # Convert 1-based idx to 0-based
    associations = np.zeros(n)
    for j in range(n):
        if j == pidx:
            associations[j] = 0.0  # Distance to self is 0
        else:
            a, b = min(pidx, j), max(pidx, j)
            pair_idx = a * (n - 1) - (a * (a - 1)) // 2 + (b - a - 1)
            associations[j] = m[pair_idx]
    return associations


# Step 7: Extract and scale associations
estimates = []
for idx in tqdm(cogtexIdxs, desc="Extracting Associations"):
    associations = getIndexPairs_numpy(idx, m_python, n)
    scaled_associations = associations / 100  # Scale as in the original script
    estimates.append(scaled_associations)
    
    
idx = cogtexIdxs[0]

#print(idx)


#associations = getIndexPairs_numpy(idx, m_python, n)
#print(f"Raw associations: {associations}")
#print(f"Scaled associations: {associations / 100}")
#print(m_python)    

# Step 8: Convert to pandas Series with gene symbols
gene_symbols = x['Gene Symbol'].tolist()
estimates_py = [pd.Series(est, index=gene_symbols) for est in estimates]


# Step 9: Create a dictionary with query genes as keys
estimates_dict_2nd = dict(zip(query_list_2ndorder, estimates_py))
    

queried_genes = set(estimates_dict.keys())
third_order_genes = set()


for gene in tqdm(second_order_genes, desc="Processing 3rd-Order Genes", leave=False):
    if gene in estimates_dict_2nd:
        est = estimates_dict_2nd[gene]
        if use_absolute:
            condition = (est.abs() > threshold) & (est.index != gene)
        else:
            condition = (est > threshold) & (est.index != gene)
        
        coexpressed_genes = est[condition].index
        third_order_candidates = set(coexpressed_genes) - queried_genes - first_order_genes - second_order_genes
        third_order_genes.update(third_order_candidates)

unique_count_3rd_order = len(third_order_genes)
print(f"Number of unique 3rd-order peripheral genes: {unique_count_3rd_order}")


####

#%%
x = [0,1,2,3]

####  r = abs(0.1)

y_01 = [3, 3997, 24731, 2470 ]
y_01_cum = [3, 4000, 28731, 31201]

###  r = abs(0.2)

### Unique
y_02 = [3, 597, 17482, 6277]
y_02_cum = [3, 600, 18082, 24359]


###  r = abs(0.30)

y_03 = [3, 61, 4241, 11202 ]
y_03_cum =[3, 64, 4305, 15507 ]

###  r = abs(0.40)

y_04 = [3, 16, 24, 4]
y_04_cum =[3, 19, 43, 47]

###  r = abs(0.50)

y_05 = [3, 8, 6, 1 ]
y_05_cum =[3, 11, 17, 18]


###  r = abs(0.60)

y_06 = [3, 3, 3,  4 ]
y_06_cum =[3, 6, 9, 13]

###  r = abs(0.70)

y_07 = [3, 1, 0, 0 ]
y_07_cum =[3, 4, 0, 0]

###  r = abs(0.80)

y_08 = [3, 0, 0, 0]
y_08_cum =[3, 0, 0 ,0 ]

###  r = abs(0.90)

y_09 = [3, 0, 0, 0]
y_09_cum =[3, 0, 0, 0]

# Plot the curves
plt.plot(x, y_01_cum, label = '|r|>0.1' )    # Curve 1
plt.plot(x, y_02_cum, label = '|r|>0.2' ) 
plt.plot(x, y_03_cum, label = '|r|>0.3' ) 
plt.plot(x, y_04_cum, label = '|r|>0.4' ) 
plt.plot(x, y_05_cum, label = '|r|>0.5' ) 
plt.plot(x, y_06_cum, label = '|r|>0.6' ) 
plt.plot(x, y_07_cum, label = '|r|>0.7' ) 
#plt.plot(x, y_08_cum, label = '|r|>0.8' ) 
#plt.plot(x, y_09_cum, label = '|r|>0.9' ) 

# Add data points
plt.scatter(x, y_01_cum, s=50, marker='o')
plt.scatter(x, y_02_cum, s=50, marker='o')
plt.scatter(x, y_03_cum, s=50, marker='o')
plt.scatter(x, y_04_cum, s=50, marker='o')
plt.scatter(x, y_05_cum, s=50, marker='o')
plt.scatter(x, y_06_cum, s=50, marker='o')
plt.scatter(x, y_07_cum, s=50, marker='o')
#plt.scatter(x, y_08_cum, s=50, marker='o')
#plt.scatter(x, y_09_cum, s=50, marker='o')

# Set y-axis to logarithmic scale
plt.yscale('log')
plt.ylim(-100, 70000)

# Customize x-axis ticks
plt.xticks(ticks=[0, 1, 2, 3], labels=["0th order", "1st order", "2nd order", "3rd order"])

# Customize the plot
plt.title('Unique Peripherals for Early Onset Cardiomyopathy', fontsize = 12)
plt.xlabel('Order of Peripheral', fontsize = 12)
plt.ylabel('Cumulative Number of Unique Peripherals', fontsize = 10.5)
plt.legend()
plt.grid(True, which="major", ls="--")  # Grid for both major and minor ticks

# Display the plot
plt.show()
#%% Imports, CSVs and dictionary setups

"""Created by Shaked on 14.03.2023
Delta distances between CA of subunits"""

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math
import re
import chain_contactMap as condif
import csv

"""Import chain .dat files"""
datA = '/home_c/shaked/ftsk/chainA/chainA.dat'
datB = '/home_c/shaked/ftsk/chainB/chainB.dat'
datC = '/home_c/shaked/ftsk/chainC/chainC.dat'
datD = '/home_c/shaked/ftsk/chainD/chainD.dat'
datE = '/home_c/shaked/ftsk/chainE/chainE.dat'
datF = '/home_c/shaked/ftsk/chainF/chainF.dat'

"""Read CSVs for chains"""
csvA = pd.read_csv('/home_c/shaked/ftsk/chainA/chainAdf.csv')
csvB = pd.read_csv('/home_c/shaked/ftsk/chainB/chainBdf.csv')
csvC = pd.read_csv('/home_c/shaked/ftsk/chainC/chainCdf.csv')
csvD = pd.read_csv('/home_c/shaked/ftsk/chainD/chainDdf.csv')
csvE = pd.read_csv('/home_c/shaked/ftsk/chainE/chainEdf.csv')
csvF = pd.read_csv('/home_c/shaked/ftsk/chainF/chainFdf.csv')
print(len(csvD['i']))

"""Dictionaries"""
dicA = condif.dicA
dicB = condif.dicB
dicC = condif.dicC
dicD = condif.dicD
dicE = condif.dicE
dicF = condif.dicF
# print(dicF)

#%% Chain distance arrays (from Traj file) and Heatmap for each chain
"""Create arrays from chain .dat file, i long and j wide (same)
for subunits - size is of each subunit (subunits used to be different).
Uses Trajectories and calculates via np.linalg.norm"""

def create_distances_array(dat):
    with open(dat, 'r') as FILEdat:
        dat = FILEdat.readlines()
        xyz = {}
        for line in dat:
            if 'size of chain' in line:
                # print(line)
                bead_no = int(re.search(r'\d+', line).group())
                distances_array = np.zeros((bead_no, bead_no))
            elif 'size of chain' not in line:
                if "CA" in line:
                    line = line.split()
                    # print(line)
                    xyz[int(line[0])] = [float(line[4]), float(line[5]), float(line[6])]
            else:
                continue
        for i in range(1, bead_no + 1):
            for j in range(i, bead_no):
                coord_i = np.array(xyz[i])
                coord_j = np.array(xyz[j])                
                distances_array[i][j] = np.linalg.norm(coord_i - coord_j)
                distances_array[j][i] = np.linalg.norm(coord_i - coord_j)
        # return xyz
        # return bead_no
                               
    return distances_array    

chainA_distances = create_distances_array(datA)
chainB_distances = create_distances_array(datB)
chainC_distances = create_distances_array(datC)
chainD_distances = create_distances_array(datD)
chainE_distances = create_distances_array(datE)
chainF_distances = create_distances_array(datF)

"""Just a check - in order to reach a certain index in the protein,
315 needs to be subtracted from the desired index"""
check = chainF_distances[384-314][670-314]
print(check)
print('Heatmap of chain distances')
sns.heatmap(chainC_distances)

"""lengths of chains in X_beads"""
A_beads = (len(chainA_distances))
B_beads = (len(chainB_distances))
C_beads = (len(chainC_distances))
D_beads = (len(chainD_distances))
E_beads = (len(chainE_distances))
F_beads = (len(chainF_distances))

"""Chains used to be different. Now they are the same so this is irrelevant"""
#%% C_F - delta differences array (from Traj file)
"""Create array of delta distance differences between chains for each pair of beads"""
def distance_differences(arr1, len1, arr2, len2):
    array_size = min(len1, len2)
    differences_array = np.zeros((array_size, array_size))
    # print('Array size: ', differences_array.shape)
    for i in range(1, array_size + 1):
        for j in range(i, array_size):
            differences_array[i][j] = abs(arr1[i][j] - arr2[i][j])
            differences_array[j][i] = abs(arr1[i][j] - arr2[i][j])
    return differences_array

C_F = distance_differences(chainF_distances, F_beads, chainC_distances, C_beads)
# print(C_F)

"""Check indices"""
# print(chainC_distances[330][335])
# print(C_F[383-314][657-314])

"""Save array to DF then CSV"""
C_F_df = pd.DataFrame(C_F)
C_F_df.to_csv(r'/home_c/shaked/ftsk/C_F_df.csv')

# A_D = distance_differences(chainA_distances, A_beads, chainD_distances, D_beads)

print('Heatmap of delta distances between chains C and F')
map = sns.heatmap(C_F, cmap='PuRd')
map.invert_yaxis()

# """Filter heatmap to contain only distances larger than x Angstrom"""
# C_F_filtered = np.where(abs(C_F) <= 11 , 0, C_F)
# print('Filtered heatmap of delta distances between chains C and F')
# map = sns.heatmap(C_F_filtered, cmap='PuRd')
# map.invert_yaxis()

#%% Scatter plot from the above filtered heatmap
print('Scatter plot of delta distances between two chains')

"""Get dimensions of array"""
nrows, ncols = C_F.shape            
rows, cols = np.meshgrid(np.arange(nrows) + 315, np.arange(ncols) + 315)            # create a meshgrid of row and column indices

"""create a boolean mask for values greater than x"""
mask = C_F >= 11
"""apply the mask to the row indices and flatten"""
x = rows[mask]
"""apply the mask to the column indices and flatten"""
y = cols[mask]
"""apply the mask to the color values and flatten"""
c = C_F[mask]

plt.scatter(x, y, c=c, cmap='RdPu', vmin= 0, s = 2)
plt.xlabel('Row index')
plt.ylabel('Column index')
plt.colorbar()
plt.show()

#%%# %% Chain-unique and shared contacts
"""Define dictionaries for unique contacts of each chain and of shared contacts"""
def unique_2_chain(dic1, dic2):
    """contacts unique to chain1 and not in chain2"""
    unique1 = {}
    for key in dic1.keys():
        if key not in dic2.keys():
            unique1[key] = dic1[key]
    return unique1

uniqueC = unique_2_chain(dicC, dicF)
uniqueF = unique_2_chain(dicF, dicC)
print(str(len(uniqueC)), 'contacts in chain C and not F')
print(str(len(uniqueF)), 'contacts in chain F and not C')

def sharedCF(dic1, dic2):
    """contacts shared between chain1 and chain2"""
    shared = {}
    for key in dic1.keys():
        if key in dic2.keys():
            shared[key] = dic1[key]
    return shared

shared = sharedCF(dicC, dicF)
print(str(len(shared)), 'shared contacts between chains C and F')

# %% New filtered scatter contact delta distance ***from chain dictionaries***
print('Delta distances between chains C and F (>=20 % change)')

"""create dic of values are deltas of contact lengths from dicC and dicF
and contact pairs in both chains have a big delta distance"""
unique_and_shared = {}
"""contacts that are in either C or F and also both chains"""
for key in dicC.keys() | dicF.keys():
    if key in dicC and key in dicF:
        if abs(dicC[key] / dicF[key]) <= 0.8:
            value = C_F[key[0] - 315][key[1] - 315]
            unique_and_shared[key] = value
        elif abs(dicC[key] / dicF[key]) >= 1.2:
            value = C_F[key[0] - 315][key[1] - 315]
            unique_and_shared[key] = value
        else:
            continue
    else:
        value = abs(C_F[key[0] - 315][key[1] - 315])
        unique_and_shared[key] = value

"""create array from dic"""
bead_1 = []
bead_2 = []
values = []
for key, value in unique_and_shared.items():            
    bead_1.append(key[0])
    bead_2.append(key[1])
    values.append(value)

"""stack lists into a 3D array"""
delta_contact_distance = np.column_stack((bead_1, bead_2, values))
# print(delta_contact_distance)

plt.scatter(delta_contact_distance[:, 0], delta_contact_distance[:, 1], c=delta_contact_distance[:, 2], s=1, cmap = 'RdPu', vmin = 0)
plt.xlim([300, 800])
plt.ylim([300, 800])
plt.colorbar()
plt.show()

"""**************************************************************************"""
"""Select rows where value is greater than x"""
"""**************************************************************************"""
print('Delta distances between chains C and F (>=20 % change \n AND delta distance above 4 Angstrom)')
mask = delta_contact_distance[:,2] >= 4
significant_distance = delta_contact_distance[mask]

fig, ax = plt.subplots()
sc = ax.scatter(significant_distance[:,0], significant_distance[:,1], c=significant_distance[:,2].astype(float), s = 1, cmap='RdPu', vmin = 0)
cbar = fig.colorbar(sc)
plt.xlim([300, 800])
plt.ylim([300, 800])
plt.show()

"""Get coordinates of each plotted point"""
x_coords = significant_distance[:,0]            
y_coords = significant_distance[:,1]

"""Show contacts with significant delta distance
*****that are in the relevant sections of the protein*****"""
print('Contacts with significant delta distance:')
contacts = []
for i in range(len(x_coords + 1)):
    contacts.append('(' + str(int(x_coords[i])) + ', ' + (str(int(y_coords[i])))+ ')')
    pairs = [eval(pair_str) for pair_str in contacts if isinstance(eval(pair_str), tuple)]
filtered_pairs = [pair for pair in pairs if int(pair[0]) <= 500 or int(pair[0]) >=600]
filtered_pairs.append('Total number of contacts: ' + str(len(filtered_pairs)))
print(filtered_pairs)

"""Show which contacts were added to delta distance list.
Input - high threshold (deltas list in larger delta),
low threshold (deltas list in larger delta)"""
high_threshold = []
low_threshold = []
added_contacts = [con for con in low_threshold if con not in high_threshold]
print('Contacts added: ', added_contacts)

#%%% Change contacts to indices from .dat file

# print(filtered_pairs)
indices_filtered_pairs = filtered_pairs.pop()
indices_filtered_pairs = [tuple(x - 314 if isinstance(x, int) else x for x in tup) for tup in filtered_pairs]
# print(indices_filtered_pairs)
# print((indices_filtered_pairs[0]))

for item in indices_filtered_pairs:
    print(item)

#%% Is contact in chain check

def con_in_chain(key, dic1, dic2):
    if key in dic1.keys():
        if key in dic2.keys():
            return('key present in both dics')
        else:
            return('key present only in dic1')
    else:
        if key in dic2.keys():
            return('key present only in dic2')
        else:
            return ('key not found')


a = con_in_chain((385, 671), dicC, dicF)
print(a)

# %%

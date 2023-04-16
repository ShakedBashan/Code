#%% Imports, CSVs and dictionary setups------------------------------------------------------------------------------------------------------------------------------------------------------------------

# Created by Shaked on 08.03.2023
# Last changed on 12.03.2023
# Shows differences between contacts of chains

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import math

# Read CSVs for chains
chainA = pd.read_csv('/home_c/shaked/ftsk/chainA/chainAdf.csv')
chainB = pd.read_csv('/home_c/shaked/ftsk/chainB/chainBdf.csv')
chainC = pd.read_csv('/home_c/shaked/ftsk/chainC/chainCdf.csv')
chainD = pd.read_csv('/home_c/shaked/ftsk/chainD/chainDdf.csv')
chainE = pd.read_csv('/home_c/shaked/ftsk/chainE/chainEdf.csv')
chainF = pd.read_csv('/home_c/shaked/ftsk/chainF/chainFdf.csv')

# Create dictionaries for chains
dicA = chainA.set_index(['i', 'j'])['Distance'].to_dict()
dicB = chainB.set_index(['i', 'j'])['Distance'].to_dict()
dicC = chainC.set_index(['i', 'j'])['Distance'].to_dict()
dicD = chainD.set_index(['i', 'j'])['Distance'].to_dict()
dicE = chainE.set_index(['i', 'j'])['Distance'].to_dict()
dicF = chainF.set_index(['i', 'j'])['Distance'].to_dict()
print('length of dicA - ' + str(len(dicA)))
print('length of dicB - ' + str(len(dicB)))
print('length of dicC - ' + str(len(dicC)))
print('length of dicD - ' + str(len(dicD)))
print('length of dicE - ' + str(len(dicE)))
print('length of dicF - ' + str(len(dicF)))
# print(dicF)
#%% Functions

def map_only (dic1, dic2):
    """Maps contacts only in chain 1 and not 2"""
    only_chain1 = {(k): [dic1[k]] for k in dic1 if k not in dic2}
    print('Unique chain contacts that appear in chain1 and not chain2: ', len(only_chain1))
    chain1_i, chain1_j = zip(*only_chain1.keys())
    plt.scatter(chain1_i, chain1_j, s=1, color = 'blue')
    plt.show()


def map_commons (dic1, dic2):
    """Maps common contacts between two chains"""
    both_chains = {(k): [dic1[k], dic2[k]] for k in dic1 if k in dic2}
    # print(both_chains)
    print('Common contacts:', len(both_chains))
    both_i, both_j = zip(*both_chains.keys())
    plt.scatter(both_i, both_j, s=1, color = 'blueviolet')
    plt.show()


def map_changed_contacts(dic1, dic2):
    """Contacts common between chains but their distance changed dramatically"""
    same_contacts = {(k): [dic1[k], dic2[k]] for k in dic1 if k in dic2}
    same_contacts_sqrt = {k: [math.sqrt(v) for v in values] for k, values \
        in same_contacts.items()}
    threshold = 0.2
    significant_differences = {k: values for k, values in same_contacts_sqrt.items() \
        if abs(values[0] - values[1]) > threshold * (values[0] + values[1]) / 2}
    changed_i, changed_j = zip(*significant_differences.keys())
    plt.scatter(changed_i, changed_j, s=1, color = 'black')
    print('Contacts in common to both chains that have changed > 20%: ', len(significant_differences))
    plt.show()
    return significant_differences

#%% Code
"""Map differences between chains:"""
map_only(dicC, dicF)

"""Map common contacts between chains:"""
map_commons(dicC, dicF)

"""Map common contacts that changed > 20%:"""
changedCF = map_changed_contacts(dicC, dicF)
print(changedCF)


"""Create maps that can be shown together"""
only_chainA = {(k): [dicA[k]] for k in dicA if k not in dicF}
# print(only_chainA)
print('Chain A unique contacts:', len(only_chainA))
Ai, Aj = zip(*only_chainA.keys())
plt.scatter(Ai, Aj, s=1, color = 'blue')
plt.show()

only_chainB = {(k): [dicB[k]] for k in dicB if k not in dicF}
# print(only_chainB)
print('Chain B unique contacts:', len(only_chainB))
Bi, Bj = zip(*only_chainB.keys())
plt.scatter(Bi, Bj, s=1, color = 'cornflowerblue')
plt.show()

only_chainC = {(k): [dicC[k]] for k in dicC if k not in dicF}
# print(only_chainC)
print('Chain C unique contacts:', len(only_chainC))
Ci, Cj = zip(*only_chainC.keys())
plt.scatter(Ci, Cj, s=1, color = 'lime')
plt.show()

only_chainD = {(k): [dicD[k]] for k in dicD if k not in dicF}
# print(only_chainD)
print('Chain D unique contacts:', len(only_chainD))
Di, Dj = zip(*only_chainD.keys())
plt.scatter(Ci, Cj, s=1, color = 'gold')
plt.show()

only_chainE = {(k): [dicE[k]] for k in dicE if k not in dicF}
print(only_chainE)
print('Chain E unique contacts:', len(only_chainE))
Ei, Ej = zip(*only_chainE.keys())
plt.scatter(Ei, Ej, s=1, color = 'darkorange')
plt.show()

only_chainF = {(k): [dicF[k]] for k in dicF if k not in dicC}
print(only_chainF)
print('Chain F unique contacts:', len(only_chainF))
Fi, Fj = zip(*only_chainF.keys())
plt.scatter(Fi, Fj, s=1, color = 'red')
plt.show()

"""Show unique chain contacts on the same plot:"""
# plt.scatter(Ai, Aj, s=1, color = 'blue')
# plt.scatter(Bi, Bj, s=1, color = 'cornflowerblue')
plt.scatter(Ci, Cj, s=1, color = 'lime')
# plt.scatter(Di, Dj, s=1, color = 'gold')
# plt.scatter(Ei, Ej, s=1, color = 'darkorange')
plt.scatter(Fi, Fj, s=1, color = 'red')
plt.show()
#%% Runs

changedCF = map_changed_contacts(dicC, dicF)
print(changedCF)

map_only(dicC, dicF)
map_only(dicF, dicC)
#%% OLD DRAFTS

# def differences_from_dictionaries(dic1, dic2):
#     same_contacts = {(k): [dic1[k], dic2[k]] for k in dic1 if k in dic2}
#     same_contacts_sqrt = {k: [math.sqrt(v) for v in values] for k, values \
#         in same_contacts.items()}
#     threshold = 0.2
#     significant_differences = {k: values for k, values in same_contacts_sqrt.items() \
#         if abs(values[0] - values[1]) > threshold * (values[0] + values[1]) / 2}
#     # for key in dic1:
#     #     if key not in contact_changes.keys():
#     #         significant_differences[key[0]] = key[1]
#     # for key in dic2:
#     #     if key not in contact_changes.keys():
#     #         significant_differences[key[0]] = key[1]
#     # return significant_differences
#     return significant_differences


# significant_differences = differences_from_dictionaries(dicC, dicF)
# print(significant_differences)
# print(len(significant_differences))


# if '386' in significant_differences:
#     print(significant_differences.key())
# else:
#     print('not found')

# x = list(significant_differences.keys())
# y = list(significant_differences.values())

# plt.scatter(x, y, s=10)

# plt.xlabel('X Axis')
# plt.ylabel('Y Axis')
# plt.show()


# def create_map(df_chain):
#     plt.scatter(df_chain['i'], df_chain['j'], marker=("."), color = 'darkred')
#     plt.xlabel('i')
#     plt.ylabel('j')
#     return plt.show()




#%% END
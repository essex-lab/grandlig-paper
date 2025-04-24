# Modules to import
import argparse
import os

import scipy.stats
from simtk.unit import *
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.optimize import curve_fit
import pickle
import re

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-f", '--folders', help='Input folder names of the fragment repeats', type=str, nargs='+')
parser.add_argument("-fn", '--name', help='Input file name of the log file', type=str)
parser.add_argument("-frag", '--frag', help='Input name of the fragment', type=str, default='Ligand')

args = parser.parse_args()


# Data needs to be in form B, repeat, n_moves, avg_Occ


def get_avgN(file):
    avgN_nmoves = []
    print(file)
    with open(file, 'r') as fi:
        for line in fi.readlines():
            if 'completed' in line:
                try:
                    avg_N = float(line.split()[-9].strip('.'))
                except ValueError:
                    avg_N = float(line.split()[-1])
                avgN_nmoves.append(avg_N)

    return avgN_nmoves


repeats = args.folders
n_repeats = len(repeats)
log_name = args.name

unique_Bs = []
for repeat in repeats:
    lig_name_glob = re.sub("([\[\]])", "[\\1]", repeat)
    folders = glob(f'{lig_name_glob}/*.*/')
    for f in folders:
        b = float(f.split('/')[-2])
        if b not in unique_Bs:
            unique_Bs.append(b)

Bs = unique_Bs.copy()

B_ids = []

for i in range(len(Bs)):
    B_ids.append(i)
Bs.sort()
print(Bs)

# Save the B values as a .npy file
#np.save(f"{args.frag}_B_values.npy", Bs)


avg_N_mols = np.zeros((n_repeats, len(B_ids)))  # Dict of B values and the current N waters at each frame

# Determine the maximum number of moves across all repeats and Bs
max_moves = 0
for repeat in repeats:
    for B in Bs:
        try:
            avg_N = get_avgN(f'{repeat}/{B}/{log_name}')
            if len(avg_N) > max_moves:
                max_moves = len(avg_N)
        except FileNotFoundError:
            print(f"Not found file: {repeat}/{B}/{log_name}")
            continue

# Initialize the array to store the data
full_data = np.full((len(Bs), max_moves, n_repeats), np.nan)


# Populate the array with the data
for repeat_id, repeat in enumerate(repeats):
    for i, B in enumerate(Bs):
        try:
            avg_N = get_avgN(f'{repeat}/{B}/{log_name}')
            for j, value in enumerate(avg_N):
                full_data[i, j, repeat_id] = value
        except FileNotFoundError:
            print(f"Not found file: {repeat}/{B}/{log_name}")
            continue

# Add the B values as an additional column to the full_data array
full_data_with_Bs = [Bs, full_data]
# Save the full_data_with_Bs array as a .pkl file
with open(f'{args.frag}_full_data.pkl', 'wb') as pkl_file:
    pickle.dump(full_data_with_Bs, pkl_file)
#np.save(f'{args.frag}_full_data.npy', full_data)

# Modules to import
import argparse
import os

import scipy.stats
from simtk.unit import *
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.optimize import curve_fit


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

                avg_N = float(line.split()[-9].strip('.'))
                #avg_N = float(line.split()[-1])
                avgN_nmoves.append(avg_N)

    return avgN_nmoves


repeats = args.folders
n_repeats = len(repeats)
log_name = args.name

# Get all the B values simulated for a particular ligand. N.B these dont have to be simulated in every repeat
unique_Bs = []
for repeat in repeats:
    folders = glob(f'{repeat}/*.*/')
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
np.save(f"{args.frag}_B_values.npy", Bs)


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

# Save the full_data array as a .npy file
print(full_data)
np.save(f'{args.frag}_full_data.npy', full_data)

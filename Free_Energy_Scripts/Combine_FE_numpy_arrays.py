#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to analyse the results of free energy calculations, accepts a list of saved numpy arrays containign the data
and the data is then analysed using pymbar to give the average free energy, standard error of the mean, overlap matrix and
the convergence plots
- Will Poole
"""

# Modules to import
import argparse
import numpy as np

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-u", "--energies", help="Input list of .npy files containing the free energy data", default=[], nargs='+')
args = parser.parse_args()


# Set up variables
files = args.energies

# Sort the files out based on lambda num
lambdas = []
for file in files:
    lam = int(file.split('/')[0])
    lambdas.append(lam)

lambdas, files = zip(*sorted(zip(lambdas, files)))
U = np.load(files[0])
n_lambdas = U.shape[0]
n_samples = U.shape[2]
if n_samples > 1000:
    sli = n_samples - 1000
    U = U[:, :, sli:]

print(U.shape)
print(f'Number of files to combine: {len(files)}. Check this is equal for the number of lambdas simulated ({n_lambdas})')

# Assume that the min and max are correct - maybe add an argument for this but hey ho
for i in range(n_lambdas):
    if i not in lambdas:
        print(f'Lambda = {i} not in list')


if len(files) != n_lambdas:
    raise Exception("Number of Files not Equal to number of Lambdas.")

for file in files[1:]:
    U_new = np.load(file)
    n_samples = U_new.shape[2]
    if n_samples > 1000:
        sli = n_samples - 1000
        U_new = U_new[:, :, sli:]

    U += U_new

np.save('Combined_Matrix.npy', U)

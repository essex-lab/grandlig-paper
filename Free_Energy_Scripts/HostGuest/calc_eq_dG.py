# Script to calculate dG for 3 repeats of a hydration fE of a molecule using alchemb.


import pickle
import os.path
from os.path import join
import pandas as pd

import argparse
from pymbar.timeseries import statisticalInefficiency, subsampleCorrelatedData, detectEquilibration
from pymbar import MBAR
import numpy as np
import matplotlib.pyplot as plt
import os
from openmm.unit import *
from tqdm import tqdm


lambdanr = 40
prod_time = 4000 # ps
repeats = [1, 2] #, 3]

WORK_DIR = os.getcwd()
ligand = WORK_DIR.split('/')[-2]
state = WORK_DIR.split('/')[-1] 
# print(WORK_DIR)
# print(ligand)
n_lambdas = 40
n_samples = 1000

lambdas = np.linspace(0, 1, n_lambdas)

kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * 298 * kelvin  # Define kT
kT_kcal = kT.in_units_of(kilocalories_per_mole)

dg_repeats = []
u_kln_list = []
overlap = np.zeros([n_lambdas, n_lambdas])
for repeat in repeats:
    os.chdir(f'repeat_{repeat}')
    U_kln = np.zeros((n_lambdas, n_lambdas, n_samples))
    U_kln = np.load("Combined_Matrix.npy")

    # Perform pymbar analysis
    N_k = np.zeros(n_lambdas, np.int32)
    for i in range(n_lambdas):
        A_t = U_kln[i, i, :]
        n_equil, g, neff_max = detectEquilibration(A_t)  # Detect equilibrated states for this lambda
        indicies = subsampleCorrelatedData(A_t, g=g)  # identify uncorrelated data
        # if len(indicies) < 500:
        #   print(i, len(indicies))
        N_k[i] = len(indicies)

        U_kln[i, :, 0:N_k[i]] = U_kln[i, :, indicies].T

    mbar = MBAR(U_kln, N_k, maximum_iterations=10000)
    #overlap += mbar.computeOverlap()["matrix"]
    [deltaG_ij, ddeltaG_ij, theta_ij] = mbar.getFreeEnergyDifferences(return_theta=True)
    for i in range(n_lambdas):
        dG_i = (deltaG_ij[0, i] * kT_kcal).in_units_of(kilocalorie_per_mole)
        # print('Free energy ({:.3f} -> {:.3f}) = {}'.format(lambdas[0], lambdas[i], dG_i))

    dg = deltaG_ij[0, -1] * kT_kcal
    dg_repeats.append(dg._value)
    os.chdir(WORK_DIR)

# overlap /= len(repeats)  # Average the overlap data
# plot_overlap(overlap, ligand)

# make_fe_trace(u_kln_list, kT_kcal, lambdas)

# print(dg_repeats)

mean_dg = np.mean(dg_repeats)
std_error = np.std(dg_repeats) / np.sqrt(len(dg_repeats))
print(f"{ligand}_{state}: {mean_dg:.3f} +/- {std_error:.3f}")

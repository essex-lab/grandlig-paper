# Modules to import
import argparse
import os
import pandas as pd
from simtk.unit import *
import numpy as np
import matplotlib.pyplot as plt
from openmmtools.constants import STANDARD_STATE_VOLUME
from tqdm import tqdm
import pickle
import functions
from functions import fit_curve, plot_mean_curve


def plot_convergance(steps, b50s, sems, ax):
    ax.errorbar(
        steps,
        b50s,
        yerr=sems,
        fmt="o",
        ecolor="r",
        markersize=10,
        markeredgewidth=3,
        markerfacecolor="white",
        markeredgecolor="red",
        linestyle="none",
        elinewidth=2,
    )
    ax.plot(steps, b50s, ls="-", label="Fit line", lw=2, c="r")
    ax.fill_between(
        steps, b50s[-1] - sems[-1], b50s[-1] + sems[-1], color="gray", alpha=0.2
    )

    # Add labels and legend
    ax.set_xlabel("Number of GCNCMC Moves")
    ax.set_ylabel("B50")

small_font = 16
medium_font = 18
large_font = 20

plt.rc('figure', titlesize=large_font)
plt.rc('font', size=small_font)
plt.rc('axes', titlesize=small_font)
plt.rc('axes', labelsize=medium_font)
plt.rc('xtick', labelsize=small_font)
plt.rc('ytick', labelsize=small_font)
plt.rc('legend', fontsize=small_font)

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument(
    "-r",
    "--sphere_rad",
    help="Input GCMC sphere radius is angstroms.",
    type=float,
    default=8.0,
)
parser.add_argument(
    "--mu_file",
    help="Input mu file to get muex from",
    type=str,
    default=None
)
args = parser.parse_args()

rad = args.sphere_rad * angstrom  
T = 298 * kelvin
kT = (BOLTZMANN_CONSTANT_kB * T * AVOGADRO_CONSTANT_NA).in_units_of(
    kilocalories_per_mole
)
beta = 1 / kT


try:
    mu_df = pd.read_csv(
        args.mu_file,
        converters={"Name": str},
    )
except:
    mu_df = pd.read_csv(
        args.mu_file, delim_whitespace=True,
        converters={"Name": str},
    )

ligands = mu_df["Name"].to_list()
# ligands = ["Toluene"]

df = pd.DataFrame()
df['Steps'] = np.arange(100, 1600, 100)
for i, lig in enumerate(ligands):
    fig = plt.figure(figsize=(20, 10))
    axes = fig.subplot_mosaic(
        """AA
            BB"""
    )  # Raw, Convergance    
    
    mu_ex = mu_df[mu_df["Name"] == lig]["Mu_Ex"].values[0] * kilocalories_per_mole
    mu_err = mu_df[mu_df["Name"] == lig]["Mu_Ex_err"].values[0]
    try:
        with open(f"{lig}_full_data.pkl", "rb") as f:
            full_data = pickle.load(f)
    except FileNotFoundError:
        print(f"File not found for {lig}.")
        continue

    Bs = np.asarray(full_data[0])
    full_data = full_data[1]  
    full_data[full_data[:, :, :] > 1] = 1  # Bring data greater than 1 to N_max = 1
    n_repeats = full_data.shape[-1]

    params_repeats = []
    data = []
    for j in range(n_repeats):
        #print(full_data[:, -1, j])
        for b in range(len(Bs)):
            try:
                last_idx_for_b = np.argwhere(~np.isnan(full_data[b, :, j]))[-1][0]
                full_data[b, -1, j] = full_data[b, last_idx_for_b, j]
            except IndexError:
                continue
        print(lig, j+1)
        params_repeats.append(fit_curve(Bs, full_data[:, -1, j]))

        data.append([Bs, full_data[:, -1, j]])
    params_repeats = np.asarray(params_repeats)
    mean_B50, sem_B50, dF_trans, dF_trans_err, kd, dG, dG_err, dG_repeats = functions.calc_fes_from_params(params_repeats, mu_ex, mu_err, rad)
    fe_data = (mean_B50, sem_B50, dF_trans, dF_trans_err, kd, dG, dG_err)
    
    plot_mean_curve(Bs, data, params_repeats, fe_data, axes["A"], "Raw Data")

    if i == 0 :
        with open("Titration_Results_ConvScript.txt", 'w') as fo:
            fo.write(f"{dG._value:.2f}\t{dG_err:.2f}\t{lig}")
    else:
        with open("Titration_Results_ConvScript.txt", 'a') as fo:
            fo.write(f"\n{dG._value:.2f}\t{dG_err:.2f}\t{lig}")

    # convergance bits
    steps = np.arange(100, 1600, 100)
    b50s = []
    sems = []
    for step in steps:
        data_ = full_data.copy()
        step_id = step -1
        params_repeats = []

        for j in range(n_repeats):
            for b in range(len(Bs)):
                try:
                    last_idx_for_b = np.argwhere(~np.isnan(data_[b, :step_id, j]))[-1][0]
                    data_[b, step_id, j] = data_[b, last_idx_for_b, j]
                except IndexError:
                    pass

            step_data = data_[:, step_id, j]
            params_repeats.append(fit_curve(Bs, step_data))
        
        params_repeats = np.asarray(params_repeats)
        mean_B50, sem_B50, dF_trans, dF_trans_err, kd, dG, dG_err, dG_repeats = functions.calc_fes_from_params(params_repeats, mu_ex, mu_err, rad)
        b50s.append(mean_B50)
        sems.append(sem_B50)
    df[f"B50_{lig}"] = b50s
    df[f"B50_sem_{lig}"] = sems
    plot_convergance(steps, b50s, sems, axes["B"])
    fig.savefig(f"{lig}_convergence_plot.pdf", bbox_inches='tight')
df.to_csv("Convergence_B50s.csv", index=False)
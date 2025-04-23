"""
Main script to take in occupancy data from a GCNCMC simulation and plot a titration curve. This script is intened to be used for multiple ligands with an indivdiual
 titration curve for each ligand plotted. Curves are plotted on the B and concentration scale. A third curve is plotted after boootstrapping the data between repeats.
A numpy array with the Adams values, mean fitted curve parameters, error in curve parameters and the free enegry data for each repeat is written out to a .npy file for future processing and plotting.
"""
# Modules to import
import argparse
import os

import pandas as pd
import scipy.stats
from simtk.unit import *
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.optimize import curve_fit
from openmmtools.constants import STANDARD_STATE_VOLUME
from tqdm import tqdm
from functions import fit_curve, plot_mean_curve
import functions

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
parser.add_argument("-n", '--n_boot', help='Input number of times to bootstrap.', type=int, default=1000)

parser.add_argument(
    "--mu_file",
    help="Input mu file to egt muex from",
    type=str,
    default=None
)
args = parser.parse_args()

nboot = args.n_boot

rad = args.sphere_rad * angstrom
T = 298 * kelvin
kT = (BOLTZMANN_CONSTANT_kB * T * AVOGADRO_CONSTANT_NA).in_units_of(kilocalories_per_mole)
beta = 1/kT

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

for i, lig in enumerate(ligands):
    # Setup figure
    fig = plt.figure(figsize=(15, 5))
    axes = fig.subplot_mosaic("""ABC""")

    lig_name = lig
    file_name = lig+"_NOpost_procress_data.csv"

    mu_ex = mu_df[mu_df["Name"] == lig_name]["Mu_Ex"].values[0] * kilocalories_per_mole
    mu_err = mu_df[mu_df["Name"] == lig]["Mu_Ex_err"].values[0]

    # Read the CSV file with the occupancy data
    df = pd.read_csv(file_name, header=None)
    num_columns = df.shape[1]
    n_repeats = num_columns - 1

    # Assign column names
    column_names = ["B"] + [f"r{i}" for i in range(1, num_columns)]
    df.columns = column_names

    # Filter out columns that are not r1, r2, r3, r4
    columns_to_keep = ["B", "r1", "r2", "r3", "r4"]
    n_repeats = 4
    df = df[columns_to_keep]

    data = []
    log_conc_data = []

    Bs_repeats = []
    log_concs_repeats = []

    occs_repeats = []
    params_repeats = []
    log_concs_params_repeats = []

    # Extract B data, occupancy data and append to list. Calculate curve fit too
    for j in range(1, n_repeats + 1):
        df[f"r{j}"] = (
            df[f"r{j}"].astype(str).str.strip().replace("nan", np.nan).astype(float)
        )
        df_r = df[~df[f"r{j}"].isna()]
        data.append([df["B"].values, df[f"r{j}"].values])
        log_conc_data.append(
            [np.log10(functions.calc_c_from_B(df["B"].values, mu_ex, rad)._value), df[f"r{j}"].values]
        )

        Bs_repeats.append(df_r[f"B"].values)
        log_concs_repeats.append(
            np.log10(functions.calc_c_from_B(df_r[f"B"].values, mu_ex, rad)._value)
        )

        occs_repeats.append(df_r[f"r{j}"].values)

        params_repeats.append(fit_curve(df_r[f"B"].values, df_r[f"r{j}"].values)) 
        log_concs_params_repeats.append(
            fit_curve(
                np.log10(functions.calc_c_from_B(df_r[f"B"].values, mu_ex, rad)._value),
                df_r[f"r{j}"].values,
            )
        )

    # Now data extracted, calculate the free energies
    params_repeats = np.asarray(params_repeats)
    mean_B50, sem_B50, dF_trans, dF_trans_err, kd, dG, dG_err, dG_repeats = functions.calc_fes_from_params(params_repeats, mu_ex, mu_err, rad)
    fe_data = (mean_B50, sem_B50, dF_trans, dF_trans_err, kd, dG, dG_err)

    Bs = np.asarray(df["B"])  # Get all the B values for all 
    log_concs = np.log10(functions.calc_c_from_B(Bs, mu_ex, rad)._value)

    # Calculate fitted curves
    mean_params = np.mean(params_repeats, axis=0)
    sem_params = np.std(params_repeats, axis=0) / np.sqrt(len(params_repeats))

    params_to_save = [Bs, mean_params, sem_params, dG_repeats]
    np.save(f"{lig}_mean_plot.npy", params_to_save)

    if i == 0 :
        with open("Titration_ResultsNEW.txt", 'w') as fo:
            fo.write(f"{dG._value:.2f}\t{dG_err:.2f}\t{lig}")
        with open("Titration_Results_B50.txt", "w") as fo:
            fo.write(f"Ligand\tMu_ex\tMu_ex_err\tB50\tB50_err\n")
            fo.write(f"{lig}\t{mu_ex._value}\t{mu_err}\t{mean_B50:.2f}\t{sem_B50:.2f}")

    else:
        with open("Titration_ResultsNEW.txt", 'a') as fo:
            fo.write(f"\n{dG._value:.2f}\t{dG_err:.2f}\t{lig}")
        with open("Titration_Results_B50.txt", "a") as fo:
            fo.write(f"\n{lig}\t{mu_ex._value}\t{mu_err}\t{mean_B50:.2f}\t{sem_B50:.2f}")

    plot_mean_curve(Bs, data, params_repeats, fe_data,  axes["A"], "Raw Data")

    # Now for conc scale
    plot_mean_curve(
        log_concs,
        log_conc_data,
        log_concs_params_repeats,
        fe_data,
        axes["B"],
        "Raw Data - Conc",
        scale="logconc",
    )
    plt.tight_layout()
    # plt.show()

    # Now bootstrapped plot
    data = np.asarray(data)
    bootoccs = functions.generate_bootstrap_data(nboot, data[:, 1, :].T, Bs)

    booted_params_repeats = np.asarray([fit_curve(Bs, bootoccs[:, j]) for j in range(nboot)])
    
    mean_B50, sem_B50, dF_trans, dF_trans_err, kd, dG, dG_err, dG_repeats = (
        functions.calc_fes_from_params(booted_params_repeats, mu_ex, mu_err, rad)
    )
    fe_data = (mean_B50, sem_B50, dF_trans, dF_trans_err, kd, dG, dG_err)
    plot_mean_curve(
        Bs, data, booted_params_repeats, fe_data, axes["C"], "Bootstrapped Data"
    )

    if i == 0:
        with open("Titration_ResultsNEW_booted.txt", "w") as fo:
            fo.write(f"{dG._value:.2f}\t{dG_err:.2f}\t{lig}")

    else:
        with open("Titration_ResultsNEW_booted.txt", "a") as fo:
            fo.write(f"\n{dG._value:.2f}\t{dG_err:.2f}\t{lig}")

    plt.savefig(
        f"{lig_name}-B_C_Boot.pdf", bbox_inches="tight", pad_inches=0.1
    )

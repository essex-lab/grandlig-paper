# Modules to import
import argparse
import os

import scipy.stats
from simtk.unit import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from glob import glob
from scipy.optimize import curve_fit
import pickle
import pandas as pd


small_font = 16
medium_font = 18
large_font = 20

plt.rc('figure', titlesize=large_font)
plt.rc('font', size=small_font)
plt.rc('axes', titlesize=medium_font)
plt.rc('axes', labelsize=medium_font)
plt.rc('xtick', labelsize=small_font)
plt.rc('ytick', labelsize=small_font)
plt.rc('legend', fontsize=14)

def sigmoid(x, x0, k):
    """
    1 / 1 + exp(-k(B-betadF))
    Parameters
    ----------
    x
    x0
    k

    Returns
    -------

    """
    y = 1 / (1 + np.exp(-k*(x-x0)))
    return (y)

T = 298 * kelvin
kT = (BOLTZMANN_CONSTANT_kB * T * AVOGADRO_CONSTANT_NA).in_units_of(kilocalories_per_mole)
beta = 1/kT

def calc_c_from_B(B, HFE, sphere_rad):
    c = np.exp(B - (beta*HFE)) / (AVOGADRO_CONSTANT_NA * calc_VGCMC(sphere_rad))
    return c.in_units_of(molar)

def calc_VGCMC(radius):
    V = 4/3 * np.pi * radius**3
    return V


# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-mu", "--mu_file", help="Input the pickle file", type=str)
parser.add_argument("-fe", '--fe_data', help='Input protein name', type=str)
parser.add_argument("-rad", '--sphere_rad', help='Input protein name', type=float)
parser.add_argument("-pro", '--pro', help='Input output name', type=str)

parser.add_argument("-out", '--out', help='Input output name', type=str)

# parser.add_argument("-f", '--folders', help='Input folder names of the fragment repeats', type=str, nargs='+')


args = parser.parse_args()

rad = args.sphere_rad * angstrom


try:
    mu_df = pd.read_csv(
        args.mu_file,
        delimiter=",",
        converters={"Name": str},
    )
except:
    mu_df = pd.read_csv(
        args.mu_file,
        delim_whitespace=True,
        converters={"Name": str},
    )

frags = mu_df["Name"].to_list()

fe_data = pd.read_csv(args.fe_data, delim_whitespace=True, header=None, names=["FE", "FE_err", "Name"], converters={"Name": str})

merged = pd.merge(mu_df, fe_data, on=["Name"])

merged.sort_values(by="FE", ascending=True, inplace=True)

sorted_frags = merged["Name"]
all_dgs = merged["FE"].values

fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(111)
plt.xlabel('Ligand Concentration / M')
plt.ylabel('Site Occupancy')

cmap = matplotlib.cm.get_cmap("tab10")
colors = cmap(np.linspace(0, 1, len(all_dgs)))

linestyles = ["-", '--', ':']
ls_index = 0
linestyle = linestyles[ls_index]
c_index = 0
for h, frag in enumerate(sorted_frags):
    lig_df = merged[merged['Name'] == frag]
    mu = lig_df["Mu_Ex"].values[0] * kilocalories_per_mole
    dG_mean = lig_df["FE"].values[0]
    dG_std = lig_df["FE_err"].values[0]
    _ = np.load(f"{frag}_mean_plot.npy", allow_pickle=True)
    
    Bs = _[0]
    B_fit = np.linspace(min(Bs)-10, max(Bs)+10, 300)
    mean_params = _[1]
    sem_params = _[2]


    N_fit_mean = sigmoid(B_fit, *mean_params)
    C_fit = calc_c_from_B(B_fit, mu, rad)
    x = [np.log10(i) for i in C_fit]
    x = C_fit

    ax1.plot(x, N_fit_mean, ls=linestyle, c=colors[c_index], label=f'{frag.title()}: '+r'${:.1f} \pm {:.1f}\ kcal\ mol^{{-1}}\ $'.format(dG_mean, dG_std))
    c_index += 1
    if c_index == 10:
        c_index =0
        ls_index += 1
        linestyle = linestyles[ls_index]


plt.xlim(10**-8, 10**1)
plt.xscale("log")
plt.minorticks_off()
ax1.axhline(0.5)
handles, labels = ax1.get_legend_handles_labels()
# dGs = [float(label.split()[1]) for label in labels]
all_dgs, handles, labels = zip(*sorted(zip(all_dgs, handles, labels), key=lambda t: t[0]))
# dGs, labels = zip(*sorted(zip(dGs, labels)))

plt.legend(handles, labels, loc='upper left', ncol=2, fancybox=True, shadow=True)
plt.ylim(0.00, 1.01)

# plt.title(r'T4L99A Ligands')
# plt.title(r'$\beta$-Cyclodextrin Ligands')
plt.title(r'{} Titrations'.format(args.pro))
plt.savefig(f"{args.out}.pdf", format='pdf', bbox_inches='tight')
plt.savefig(f"/home/will/scripts_for_paper/{args.out}.pdf", format='pdf', bbox_inches='tight')

plt.savefig(f"{args.out}.png", bbox_inches='tight')

plt.show()

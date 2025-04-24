"""
Script which reads in data from titration analysis (.npy) and plots the dtitration curves for multiple ligands aswell as a bar plot detailing the calculated free energy.
"""

# Modules to import
import argparse
from simtk.unit import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from glob import glob
from scipy.optimize import curve_fit
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
plt.rc('legend', fontsize=17)

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
parser.add_argument("-mu", "--mu_file", help="Input Mu ex file", type=str)
parser.add_argument("-fe", '--fe_data', help='Input Free Energy Data', type=str)
parser.add_argument("-rad", '--sphere_rad', help='Input GCMC Sphere radius', type=float)
parser.add_argument("-pro", '--pro', help='Input protein name', type=str)
parser.add_argument("-out", '--out', help='Input output name', type=str)


# parser.add_argument("-f", '--folders', help='Input folder names of the fragment repeats', type=str, nargs='+')


args = parser.parse_args()

rad = args.sphere_rad * angstrom

# Reads in the file containing ligand names and the excess chemical potential.
try:
    mu_df = pd.read_csv(
        args.mu_file,
        delimiter=',',
        converters={"Name": str},
    )
except:
    mu_df = pd.read_csv(
        args.mu_file, delim_whitespace=True,
        converters={"Name": str},
    )


frags = mu_df["Name"].to_list()  # Gets the names of the fraagments

fe_data = pd.read_csv(args.fe_data, delim_whitespace=True, header=None, names=["FE", "FE_err", "Name"], converters={"Name": str})  # Get the final averaged free energy estimates from a file

merged = pd.merge(mu_df, fe_data, on=["Name"]) # Merge the two dataframes on the name column

merged.sort_values(by="FE", ascending=True, inplace=True)  # Sort the dataframe by the free energy estimates

sorted_frags = merged["Name"]
all_dgs = merged["FE"].values
all_sems = merged["FE_err"].values
all_dGs_repeats = []

sorted_dGs = all_dgs

fig = plt.figure(figsize=(20, 10))
axes = fig.subplot_mosaic("""AB""")
axes["A"].set_xlabel('Ligand Concentration / M')
axes["A"].set_ylabel('Site Occupancy')

cmap = matplotlib.cm.get_cmap("tab20")
cmap2 = matplotlib.cm.get_cmap("tab20b")

colors = cmap.colors + cmap2.colors

linestyles = ["-", '--', ':']
ls_index = 0
linestyle = linestyles[ls_index]
c_index = 0
fitted_data = pd.DataFrame()
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
    dg_repeats = _[3]
    all_dGs_repeats.append(dg_repeats)
    N_fit_mean = sigmoid(B_fit, *mean_params)

    fitted_data[f"Bs_{frag}"] = B_fit
    fitted_data[f"Mean_fit_{frag}"] = N_fit_mean

    C_fit = calc_c_from_B(B_fit, mu, rad)
    x = [np.log10(i) for i in C_fit]
    x = C_fit

    axes["A"].plot(x, N_fit_mean, ls=linestyle, c=colors[c_index], label=f'{frag.title()}') #: ' + r'$\tau={:.2f}$'.format(tau))
    c_index += 1
    if c_index == 10:
       c_index =0
       ls_index += 1
       linestyle = linestyles[ls_index]

print(all_dGs_repeats)
fitted_data.to_csv(f"{args.out}_fitted_data.csv", index=False)
axes["A"].set_xlim(10**-6, 10**1)
axes["A"].set_xscale("log")
axes["A"].axhline(0.5, c='k')
handles, labels = axes["A"].get_legend_handles_labels()
# dGs = [float(label.split()[1]) for label in labels]
all_dgs, handles, labels, all_dGs_repeats = zip(*sorted(zip(all_dgs, handles, labels, all_dGs_repeats), key=lambda t: t[0]))
# dGs, labels = zip(*sorted(zip(dGs, labels)))

# plt.legend(handles, labels, loc='upper left', ncol=1, fancybox=True, shadow=True)
axes["A"].set_ylim(0.00, 1.01)

all_dgs_value =  all_dgs
all_dgsstd_values = all_sems


axes["B"].bar(labels, all_dgs_value, yerr=all_dgsstd_values, color=colors, edgecolor='black', width=0.6, capsize=5)

for i, (label, dg_repeats) in enumerate(zip(labels, all_dGs_repeats)):
    x_positions = np.full(len(dg_repeats), i)  # Create x positions for the dots
    axes["B"].scatter(x_positions, dg_repeats, color='black', zorder=5)  # Plot the dots


box = axes["B"].get_position()
axes["B"].set_position([box.x0, box.y0 + (box.height - (box.height*0.7)), box.width, box.height*0.7])
axes["B"].set_xticks(axes["B"].get_xticks(), labels, rotation=90, fontsize=14)

axes["B"].set_ylabel('Free Energy / '+ r'$kcal mol^{-1}$')

# Shrink current axis by 20%
# box = axes["A"].get_position()
# axes["A"].set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
# axes["A"].legend(handles, labels, loc='center left', bbox_to_anchor=(1, 0.5))


# plt.title(r'T4L99A Ligands')
# plt.title(r'$\beta$-Cyclodextrin Ligands')
axes["A"].minorticks_off()

# fig.suptitle(r'{} Titrations'.format(args.pro))
plt.savefig(f"{args.out}.pdf", format='pdf', bbox_inches='tight')
# plt.savefig(f"{args.out}.png", bbox_inches='tight')

plt.show()

# Modules to import
import argparse
import os
import warnings

import pandas as pd
import scipy.stats
from simtk.unit import *
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.optimize import curve_fit, OptimizeWarning
from openmmtools.constants import STANDARD_STATE_VOLUME
import pickle

small_font = 16
medium_font = 18
large_font = 20

plt.rc("figure", titlesize=large_font)
plt.rc("font", size=small_font)
plt.rc("axes", titlesize=small_font)
plt.rc("axes", labelsize=medium_font)
plt.rc("xtick", labelsize=small_font)
plt.rc("ytick", labelsize=small_font)
plt.rc("legend", fontsize=small_font)

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-G", "--HFE", help="Input the mu_ex file")
parser.add_argument(
    "-r",
    "--sphere_rad",
    help="Input GCMC sphere radius is angstroms.",
    type=float,
    default=8.0,
)
parser.add_argument(
    "-err",
    "--plot_error",
    help="Do you want to plot the errors?.",
    default=False,
    action="store_true",
)
parser.add_argument(
    "-pro", "--protein", help="What protein is this for?", default=None, type=str
)
parser.add_argument(
    "--square",
    help="do you want the plot to be 10x10 square?",
    default=False,
    action="store_true",
)


# parser.add_argument("-f", '--folders', help='Input folder names of the fragment repeats', type=str, nargs='+')


args = parser.parse_args()

warnings.simplefilter("error", OptimizeWarning)


def calc_VGCMC(radius):
    V = 4 / 3 * np.pi * radius**3
    return V


def calc_LigVol(conc):
    V_l = 1 / (AVOGADRO_CONSTANT_NA * conc)
    return V_l.in_units_of(angstroms**3)


def calc_c_from_B(B, HFE, sphere_rad):
    c = np.exp(B - (beta * HFE)) / (AVOGADRO_CONSTANT_NA * calc_VGCMC(sphere_rad))
    return c.in_units_of(molar)


def calcB(HFE, sphere_rad, V_L):
    test = calc_VGCMC(sphere_rad) / (V_L)
    B = (beta * HFE) + np.log(test)
    return B


def sigmoid(x, x0, k):
    y = 1 / (1 + np.exp(-k * (x - x0)))
    return y


def inverse_sigmoid(Y, x0, k):
    return (-np.log((1 - Y) / Y) / k) + x0


def bootstrap(n, occupancies, Bs, mu_ex):
    booted_params = []
    n_skipped = 0
    for i in range(n):
        B_sample, occ_sample = [], []
        for j, B in enumerate(Bs):
            B_sample.append(B)
            rand_data_set = np.random.randint(0, len(occupancies))
            while np.isnan(occupancies[rand_data_set][j]):
                rand_data_set = np.random.randint(0, len(occupancies))
            # print(rand_data_set)
            occ_sample.append(occupancies[rand_data_set][j])
        # print(occ_sample)
        concs = calc_c_from_B(Bs, mu_ex, args.sphere_rad * angstroms)
        log_concs = np.log10(concs)
        # print(concs, occupancies)
        try:  # Add in initial fits and get x number of fits per bootstrap and take the best?
            initial_guess = [np.median(log_concs), 1]
            popt, pcov = curve_fit(
                sigmoid,
                log_concs,
                occ_sample,
                p0=initial_guess,
                maxfev=10000,
                method="trf",
            )
            booted_params.append(popt)
        except:
            n_skipped += 1
            continue
    #        print(f'Skipped : {n_skipped}')

    return booted_params


def calc_B_from_c(c, HFE, sphere_rad):
    B = np.log(c * (AVOGADRO_CONSTANT_NA * calc_VGCMC(sphere_rad))) + (beta * HFE)
    # print(B)
    return B


if args.square:
    fig = plt.figure(figsize=(18.5, 6))
else:
    fig = plt.figure(figsize=(20, 10))

ax1 = fig.add_subplot(111)
plt.xlabel("log(concentration)")
plt.ylabel("Site Occupancy")

colors = {
    "Crimson": "d7263d",
    "Brandy": "914542",
    "Portland Orange": "f46036",
    "Space Cadet": "2e294e",
    "Persian Green": "1b998b",
    "Malachite": "04e762",
    "June Bud": "c5d86d",
}
colors = list(colors.values())

frags = []
mu_exs = []
T = 298 * kelvin
kT = BOLTZMANN_CONSTANT_kB * T * AVOGADRO_CONSTANT_NA
beta = 1 / kT
plot_data = {}

mu_df = pd.read_csv(args.HFE)
frags = list(mu_df["Name"])
mu_exs = list(mu_df["Mu_Ex"])

all_dgs = []
for h, frag in enumerate(frags):
    mu_ex = mu_exs[h] * kilocalories_per_mole
    B_eq = (beta * mu_ex) + np.log(
        calc_VGCMC(args.sphere_rad * angstroms) / STANDARD_STATE_VOLUME
    )  # Bstd
    B_eq_kt = (B_eq * kT).in_units_of(kilocalories_per_mole)
    lig_name = frag
    file_name = f"{frag}_NOpost_procress_data.csv"
    print(lig_name)
    occs = []
    B_values = []
    with open(file_name, "r") as fi:
        lines = fi.readlines()
        for line in lines:
            split = line.split(",")
            n_repeats = len(split) - 1
            B_values.append(float(split[0]))

    with open(file_name, "r") as fi:
        lines = fi.readlines()
        for repeat in range(n_repeats):
            repeat_occs = []
            for line in lines:
                split = line.split(",")
                repeat_occs.append(float(split[repeat + 1]))
            repeat_occs = np.asarray(repeat_occs)
            occs.append(repeat_occs)

    B_values = np.asarray(B_values)

    Bs = B_values
    all_Bs = np.concatenate((Bs, np.tile(Bs, len(occs) - 1)))

    all_occs = []
    for repeat in occs:
        for i in repeat:
            all_occs.append(i)
    all_concs = calc_c_from_B(all_Bs, mu_ex, args.sphere_rad * angstroms)
    log_all_concs = np.log10(all_concs)
    tau = scipy.stats.kendalltau(log_all_concs, all_occs, nan_policy="omit")[0]

    concs = calc_c_from_B(Bs, mu_ex, args.sphere_rad * angstroms)
    log_concs = np.log10(concs)
    # Do the bootstrapping

    params = bootstrap(10000, occs, Bs, mu_ex)  # 10000
    params = np.array(params)

    N_fit_mean = np.mean(params, axis=0)
    N_fit_std = np.std(params, axis=0)

    log_kD = N_fit_mean[0]  # log(kD) by other method (better method)
    mean_Kd = (10**log_kD) * molar
    std_log_Kd = N_fit_std[0]
    std_Kd = (10**std_log_Kd) * molar

    B50 = calc_B_from_c(
        mean_Kd, mu_ex, args.sphere_rad * angstroms
    )  # also beta*df_trans
    dF_trans = (kT * B50).in_units_of(kilocalories_per_mole)

    dG_mean = (kT * np.log(mean_Kd)).in_units_of(kilocalories_per_mole)
    dG_std = (kT * np.log(std_Kd)).in_units_of(kilocalories_per_mole)

    ###### Loop to get stds
    x = np.linspace(-10, 1, 100)

    N_fit = []
    Kds = []
    for set in params:
        N_fit.append(sigmoid(x, *set))
        log_conc_50 = inverse_sigmoid(0.5, *set)
        conc_50 = 10**log_conc_50 * molar
        conc_50_micro = conc_50.in_units_of(millimolar)
        Kds.append(conc_50_micro._value)

    std_Kds = np.std(Kds)

    y = sigmoid(x, *N_fit_mean)
    y_plus = sigmoid(x, *N_fit_mean) + np.std(N_fit, axis=0)
    y_minus = sigmoid(x, *N_fit_mean) - np.std(N_fit, axis=0)

    plot_data[lig_name] = [x, y, [dG_mean, dG_std, tau]]
    # ax1.plot(x, np.quantile(N_fit, 0.5, axis=0), '-', label='{}: {:.2f} +/- {:.2f} kcal / mol '.format(lig_name, mean_dG, std_dGs)+r'($\tau={:.2f}$)'.format(tau))

    # ax1.fill_between(x, np.quantile(N_fit, 0.16, axis=0), np.quantile(N_fit, 0.84, axis=0), alpha=0.5, lw=0, color='#000080')
    # ax1.fill_between(x, np.quantile(N_fit, 0.025, axis=0), np.quantile(N_fit, 0.975, axis=0), alpha=0.1, lw=0, color='#000080')

    ax1.plot(
        x,
        y,
        "-",
        label=f"{lig_name}: "
        + r"${:.1f} \pm {:.1f}\ kcal\ mol^{{-1}}\ $".format(
            dG_mean._value, dG_std._value
        )
        + r"($\tau={:.2f}$)".format(tau),
    )
    with open("TitrationResults.txt", "a") as fo:
        fo.write("{:.1f}\t {:.1f}\t {}\n".format(dG_mean._value, dG_std._value, frag))
    all_dgs.append(dG_mean._value)
    if args.plot_error:
        ax1.fill_between(x, y_minus, y_plus, alpha=0.5, lw=0, color="#000080")
        ax1.plot(
            log_all_concs, all_occs, marker="x", linestyle="None"
        )  # , c='#'+colors[h])

# save dictionary to person_data.pkl file
with open("titration_plot_data.pkl", "wb") as fp:
    pickle.dump(plot_data, fp)
    print("dictionary saved successfully to file")

ax1.axhline(0.5)
handles, labels = ax1.get_legend_handles_labels()
# dGs = [float(label.split()[1]) for label in labels]
all_dgs, handles, labels = zip(
    *sorted(zip(all_dgs, handles, labels), key=lambda t: t[0])
)
# dGs, labels = zip(*sorted(zip(dGs, labels)))

plt.legend(
    handles, labels, loc="lower right", ncol=2, fancybox=True, shadow=True, fontsize=10
)
plt.ylim(0.00, 1.01)

plt.title(r"{} Titrations".format(args.protein))

plt.savefig(
    f"ALLFRAGSNOPOST-{args.protein}.png", bbox_inches="tight", pad_inches=0.1, dpi=600
)
plt.savefig(f"ALLFRAGSNOPOST-{args.protein}.pdf", format="pdf", bbox_inches="tight")

plt.show()

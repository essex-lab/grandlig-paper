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
    y = 1 / (1 + np.exp(-k * (x - x0)))
    return y


rad = 8 * angstrom
T = 298 * kelvin
kT = (BOLTZMANN_CONSTANT_kB * T * AVOGADRO_CONSTANT_NA).in_units_of(
    kilocalories_per_mole
)
beta = 1 / kT

def gen_text(mean_B50, sem_B50, dF_trans, dF_trans_err, dG, dG_err):
    return r"$B_{{50}} = \beta \Delta F_{{trans}} = {:.1f} \pm {:.1f}$".format(
        mean_B50, sem_B50
    ) + "\n" + r"$\Delta F_{{trans}} = {:.1f} \pm {:.1f}\ kcal\ mol^{{-1}}$".format(
        dF_trans._value, dF_trans_err._value
    ) + "\n" r"$\Delta G^{{o}}_{{bind}} = \Delta F_{{trans}} + \Delta F_{{ideal}} - \Delta G_{{sol}}$" + "\n" + r"$\Delta G^{{o}}_{{bind}} = {:.1f} \pm {:.1f}\ kcal\ mol^{{-1}}$".format(
        dG._value, dG_err
    )


def fit_curve(Bs, occs, mu_ex, rad, beta=beta):
    params = []
    Nbs, Nrs = occs.shape
    for r in range(Nrs):
        initial_guess = [np.median(Bs), 1]
        occ_list = occs[:, r]
        Bs_filtered = Bs
        occ_list_filtered = occ_list
        popt, pcov = curve_fit(
            sigmoid,
            Bs_filtered,
            occ_list_filtered,
            p0=initial_guess,
            maxfev=10000,
            nan_policy="omit",
        )
        params.append(popt)

    mean_B50 = np.mean(params, axis=0)[0]
    sem_B50 = np.std(params, axis=0)[0] / np.sqrt(len(params))

    dF_trans = (mean_B50 / beta).in_units_of(kilocalories_per_mole)
    dF_trans_err = (sem_B50 / beta).in_units_of(kilocalories_per_mole)

    dG = kT * np.log(calc_c_from_B(mean_B50, mu_ex, rad))
    return params, mean_B50, sem_B50, dF_trans, dF_trans_err, dG


def generate_bootstrap_data(n, occupancies, Bs):
    # Generate lots of randome occupancy samples

    Nbs, Nrs = occupancies.shape
    bootstrap_data = np.zeros((Nbs, n))
    for i in tqdm(range(Nbs)):
        for j in range(n):
            rand_repeat_id = np.random.randint(0, Nrs)
            # while np.isnan(occupancies[i][rand_repeat_id]):
            #     rand_repeat_id = np.random.randint(0, Nrs)
            bootstrap_data[i][j] = occupancies[i][rand_repeat_id]

    return bootstrap_data


def plot_mean_curve(Bs, occs, ax, title):
    B_fit = np.linspace(min(Bs), max(Bs), 500)
    N_fit_all = []
    tau = scipy.stats.kendalltau(Bs, np.nanmean(occs, axis=1), nan_policy="omit")[0]
    params, mean_B50, sem_B50, dF_trans, dF_trans_err, dG = fit_curve(
        Bs, occs, mu_ex, rad
    )
    dG_err = np.sqrt(dF_trans_err**2 + mu_err**2)
    # Calculate non fitted curves
    fit = []
    for p in params:
        fit.append(sigmoid(Bs, *p))

    # Calculate fitted curves
    mean_params = np.mean(params, axis=0)
    sem_params = np.std(params, axis=0) / np.sqrt(len(params))
    for p in params:
        N_fit = sigmoid(B_fit, *p)
        N_fit_all.append(N_fit)

    # Plot raw data
    Nbs, Nrs = occs.shape
    plotted_legend = False
    for r in range(Nrs):
        occ_list = occs[:, r]
        if not plotted_legend:
            ax.plot(
                Bs,
                occ_list,
                marker="x",
                linestyle="None",
                c="black",
                label="Raw data " + r"($\tau={:.3f}$)".format(tau),
            )
            plotted_legend = True
        else:
            ax.plot(Bs, occ_list, marker="x", linestyle="None", c="black")

    N_fit_mean = sigmoid(Bs, *mean_params)
    # np.std(fit, axis=0)[0] / np.sqrt(len(fit))

    B_fit = Bs
    tau_fit = scipy.stats.kendalltau(Bs, N_fit_mean, nan_policy="omit")[0]
    ax.plot(
        B_fit,
        N_fit_mean,
        "-",
        color="#04e762",
        label="Mean Fit" + r" ($\tau={:.3f}$)".format(tau_fit),
    )
    ax.fill_between(
        B_fit,
        sigmoid(Bs, *(mean_params - sem_params)),
        sigmoid(Bs, *(mean_params + sem_params)),
        alpha=0.5,
        lw=0,
        color="#000080",
    )

    # x50 = [mean_B50, mean_B50]
    # y50 = [0, 0.5]
    y_vline = sigmoid(mean_B50, *np.mean(params, axis=0))
    # print("y_vline", y_vline)

    ax.vlines(mean_B50, 0, y_vline, linestyle="--", color="k")
    # ax.plot(x50, y50, linestyle="--", c="k")

    # Add the text at the bottom right
    ax.text(
        0.98,
        0.05,
        gen_text(mean_B50, sem_B50, dF_trans, dF_trans_err, dG, dG_err),
        transform=ax.transAxes,
        verticalalignment="bottom",
        horizontalalignment="right",
        fontsize=12,
        bbox=dict(facecolor="white", alpha=0.8),
    )
    ax.set_title(title)
    ax.set_xlabel("Adams value, $B$")
    ax.set_ylabel(r"Average N")
    ax.legend(loc="upper left")


def plot_convergance(steps, b50s, sems, ax):
    ax.errorbar(
        steps,
        b50s[1:],
        yerr=sems[1:],
        fmt="o",
        ecolor="r",
        markersize=10,
        markeredgewidth=3,
        markerfacecolor="white",
        markeredgecolor="red",
        linestyle="none",
        elinewidth=2,
    )
    # ax.plot(steps, b50s[1:], marker="o", linestyle="None", mfc="none", ms=10, c='r')
    # Plot the line connecting the data points
    ax.plot(steps, b50s[1:], ls="-", label="Fit line", lw=2, c="r")
    ax.fill_between(
        steps, b50s[-1] - sems[-1], b50s[-1] + sems[-1], color="gray", alpha=0.2
    )

    # Add labels and legend
    ax.set_xlabel("Number of GCNCMC Moves")
    ax.set_ylabel("B50")


mu_df = pd.read_csv(
    args.mu_file,
    converters={"Name": str},
)

ligands = mu_df["Name"].to_list()
# ligands = ["Toluene"]
for i, lig in enumerate(ligands):
    mu_ex = mu_df[mu_df["Name"] == lig]["Mu_Ex"].values[0] * kilocalories_per_mole
    mu_err = mu_df[mu_df["Name"] == lig]["Mu_Ex_err"].values[0] * kilocalories_per_mole
    Bs = np.load(
        f"/home/will/data_6/LIG_GCNCMC_PAPER/REVIEW_ANALYSIS/TIT_CONVER/npys/{lig}_B_values.npy"
    )
    full_data = np.load(
        f"/home/will/data_6/LIG_GCNCMC_PAPER/REVIEW_ANALYSIS/TIT_CONVER/npys/{lig}_full_data.npy"
    )

    full_data[full_data[:, :, :] > 1] = 1  # Bring data to N_max = 1
    occs = full_data[:, -1, :]
    params, mean_B50, sem_B50, dF_trans, dF_trans_err, dG = fit_curve(
        Bs, occs, mu_ex, rad
    )

    dG_err = np.sqrt(dF_trans_err**2 + mu_err**2)
    print(sem_B50, mu_err, dG_err)
    if i == 0 :
        with open("Titration_Results_ConvScript.txt", 'w') as fo:
            fo.write(f"{dG._value:.2f}\t{dG_err:.2f}\t{lig}")
    else:
        with open("Titration_Results_ConvScript.txt", 'a') as fo:
            fo.write(f"\n{dG._value:.2f}\t{dG_err:.2f}\t{lig}")

    bootoccs = generate_bootstrap_data(5000, occs, Bs)

    b50s = [0]
    sems = [0]

    # convergance bits
    steps = np.arange(100, 1600, 100)
    for step in steps:
        step_data = full_data[:, step - 1, :]
        params, mean_B50, sem_B50, dF_trans, dF_trans_err, dG = fit_curve(
            Bs, step_data, mu_ex, rad
        )
        b50s.append(mean_B50)
        sems.append(sem_B50)

    fig = plt.figure(figsize=(20, 10))
    axes = fig.subplot_mosaic(
        """AB
                                CC"""
    )  # Raw, Booted, Convergance

    plot_mean_curve(Bs, occs, axes["A"], "Raw Data")
    bootoccs = generate_bootstrap_data(5000, occs, Bs)
    plot_mean_curve(Bs, bootoccs, axes["B"], "Bootstrapped")
    plot_convergance(steps, b50s, sems, axes["C"])
    fig.tight_layout()
    fig.savefig(f"/home/will/data_6/LIG_GCNCMC_PAPER/REVIEW_ANALYSIS/TIT_CONVER/plots/{lig}_convergence_plot.pdf", bbox_inches='tight')

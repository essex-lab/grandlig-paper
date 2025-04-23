"""
Script containing all functions useful for the analysis of GCNCMC simulations / titrations
- Will Poole
"""

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

# Plotting bits
small_font, medium_font, large_font = (16, 18, 20)
plt.rc("figure", titlesize=large_font)
plt.rc("font", size=small_font)
plt.rc("axes", titlesize=small_font)
plt.rc("axes", labelsize=medium_font)
plt.rc("xtick", labelsize=small_font)
plt.rc("ytick", labelsize=small_font)
plt.rc("legend", fontsize=small_font)

T = 298 * kelvin
kT = (BOLTZMANN_CONSTANT_kB * T * AVOGADRO_CONSTANT_NA).in_units_of(
    kilocalories_per_mole
)
beta = 1 / kT

# Dealing with titrations first
def calc_VGCMC(radius):
    """Calculates the volume of the GCMC sphere for a given radius

    Parameters
    ----------
    radius : float or openmm.unit.Quantity
        The radius of the sphere. If a Quantity is provided, it should have units compatible with length.

    Returns
    -------
    float or openmm.unit.Quantity
        The volume of the sphere. If the input radius is a Quantity, the output will be a Quantity with units of volume.
    """
    V = 4 / 3 * np.pi * radius**3
    return V


def calc_LigVol(conc):
    """Calculates the volume of a ligand given its concentration

    Parameters
    ----------
    conc : float or openmm.unit.Quantity
        The concentration of the ligand. If a Quantity is provided, it should have units compatible with molarity.

    Returns
    -------
    openmm.unit.Quantity
        The volume of the ligand in cubic angstroms.
    """
    V_l = 1 / (AVOGADRO_CONSTANT_NA * conc)
    return V_l.in_units_of(angstroms**3)


def calc_c_from_B(B, HFE, sphere_rad):
    c = np.exp(B - (beta * HFE)) / (AVOGADRO_CONSTANT_NA * calc_VGCMC(sphere_rad))
    return c.in_units_of(molar)


def calcB(HFE, sphere_rad, V_L):
    std = calc_VGCMC(sphere_rad) / (V_L)
    B = (beta * HFE) + np.log(std)
    return B


def sigmoid(x, x0, k):
    """
    1 / 1 + exp(-k(B-betadF))
    or
    1 / 1 + exp(-k(log10(c) - log10(Kd)))
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


def inverse_sigmoid(Y, x0, k):
    return (-np.log((1 - Y) / Y) / k) + x0


def gen_text(mean_B50, sem_B50, dF_trans, dF_trans_err, dG, dG_err):
    return r"$B_{{50}} = \beta \Delta F_{{trans}} = {:.1f} \pm {:.1f}$".format(
        mean_B50, sem_B50
    ) + "\n" + r"$\Delta F_{{trans}} = {:.1f} \pm {:.1f}\ kcal\ mol^{{-1}}$".format(
        dF_trans._value, dF_trans_err._value
    ) + "\n" r"$\Delta G^{{o}}_{{bind}} = \Delta F_{{trans}} + \Delta F_{{ideal}} - \Delta G_{{sol}}$" + "\n" + r"$\Delta G^{{o}}_{{bind}} = {:.1f} \pm {:.1f}\ kcal\ mol^{{-1}}$".format(
        dG._value, dG_err
    )


def fit_curve(Bs, occ):
    params = []
    initial_guess = [np.median(Bs), 1]

    mask = ~np.isnan(occ)
    Bs_filtered = Bs[mask]
    occ_list_filtered = occ[mask]
    if len(occ_list_filtered) < 5:
        print(f"Number of occupancy measurements is less than 5 ({len(occ_list_filtered)}). Continuing..")
        return None

    params, pcov = curve_fit(
        sigmoid,
        Bs_filtered,
        occ_list_filtered,
        p0=initial_guess,
        maxfev=50000,
        nan_policy="omit",
    )
    # params.append(popt)


    return params

def calc_fes_from_params(params, mu_ex, mu_err, rad):
    mean_B50 = np.mean(params, axis=0)[0]
    sem_B50 = np.std(params, axis=0)[0] / np.sqrt(len(params))
    # print(f"Params shape: {params.shape}")
    # print(params)
    # print(params[:, 0])

    dF_trans = (mean_B50 / beta).in_units_of(kilocalories_per_mole)
    dF_trans_err = (sem_B50 / beta).in_units_of(kilocalories_per_mole)

    kd = calc_c_from_B(mean_B50, mu_ex, rad)
    mean_dG = kT * np.log(kd)
    dG_err = np.sqrt(dF_trans_err._value**2 + mu_err**2)
    dG_repeats = kT * np.log(calc_c_from_B(params[:, 0], mu_ex, rad))

    return mean_B50, sem_B50, dF_trans, dF_trans_err, kd, mean_dG, dG_err, dG_repeats


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


def plot_mean_curve(Bs, data, params, fe_data, ax, title, scale='B'):
    mean_B50, sem_B50, dF_trans, dF_trans_err, kd, dG, dG_err = fe_data
    B_fit = np.linspace(min(Bs), max(Bs), 500)
    N_fit_all = []
    data = np.asarray(data)

    Bs_plot = data[:, 0, :]
    occs = data[:, 1, :]
    tau = scipy.stats.kendalltau(Bs, np.nanmean(occs, axis=0), nan_policy="omit")[0]

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
        Bs_list = Bs_plot[:, r]
        if not plotted_legend:
            ax.plot(
                Bs_list,
                occ_list,
                marker="x",
                linestyle="None",
                c="black",
                label="Raw data " + r"($\tau={:.3f}$)".format(tau),
            )
            plotted_legend = True
        else:
            ax.plot(Bs_list, occ_list, marker="x", linestyle="None", c="black")

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

    if scale == 'logconc':
        ax.set_xlabel(r"$\log_{10}(\mathrm{Concentration})$")

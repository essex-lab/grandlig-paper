import pandas as pd
import argparse
from pymbar.timeseries import statistical_inefficiency, subsample_correlated_data, detect_equilibration
from pymbar import MBAR
import numpy as np
import os
from openmm.unit import *
import sys
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))
from Plot_FE_convergance import make_fe_trace

def calc_mu_ex(id, n_repeats=3):
    WORK_DIR = os.getcwd()
    ligand = id
    # print(WORK_DIR)
    # print(ligand)
    repeats = range(1, n_repeats+1)
    # Load first np array to get info from

    u_ = np.load(f"{ligand}/repeat_1/lambda_0/U_matrix.npy")
    n_lambdas, _, n_samples = u_.shape
    lambdas = np.linspace(0, 1, n_lambdas)

    kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * 298 * kelvin  # Define kT
    kT_kcal = kT.in_units_of(kilocalories_per_mole)

    dg_repeats = []
    u_kln_list = []
    for repeat in repeats:
        if not os.path.isfile(f"{ligand}/repeat_{repeat}/lambda_0/U_matrix.npy"):
            continue
        U_kln = np.zeros((n_lambdas, n_lambdas, n_samples))
        for i in range(0, n_lambdas):  # Load in all the data to the big U_kln array
            npy_array = f'{ligand}/repeat_{repeat}/lambda_{i}/U_matrix.npy'
            if not os.path.isfile(npy_array):
                raise Exception(f"Cannot find file: {npy_array}")
            u_kln_ = np.load(npy_array)
            U_kln += u_kln_
        np.save(f"{ligand}/repeat_{repeat}/Combined_Ukln.npy", U_kln)

        u_kln_list.append(U_kln)

        # Perform pymbar analysis
        N_k = np.zeros(n_lambdas, np.int32)
        for i in range(n_lambdas):
            A_t = U_kln[i, i, :]
            n_equil, g, neff_max = detect_equilibration(A_t, nskip=1)  # Detect equilibrated states for this lambda
            indicies = subsample_correlated_data(A_t, g=g)  # identify uncorrelated data
            if len(indicies) < 500:
                print(i, len(indicies))
            N_k[i] = len(indicies)

            U_kln[i, :, 0:N_k[i]] = U_kln[i, :, indicies].T

        mbar = MBAR(U_kln, N_k, maximum_iterations=100000)
        [deltaG_ij, ddeltaG_ij, theta_ij] = mbar.compute_free_energy_differences(return_theta=True).values()

        for i in range(n_lambdas):
            dG_i = (deltaG_ij[0, i] * kT_kcal).in_units_of(kilocalorie_per_mole)
            print('Free energy ({:.3f} -> {:.3f}) = {}'.format(lambdas[0], lambdas[i], dG_i))

        dg = deltaG_ij[0, -1] * kT_kcal
        dg_repeats.append(dg._value)
        os.chdir(WORK_DIR)

    _, _, _ = make_fe_trace(
        u_kln_list, f"{ligand} Restraints On", kT_kcal, lambdas, f"{ligand}/{ligand}_rest_convergence.pdf"
    )
    mean_dg = np.mean(dg_repeats)
    std_error = np.std(dg_repeats) / np.sqrt(len(dg_repeats))
    print(f"{ligand}: {mean_dg:.3f} +/- {std_error:.3f}")
    return mean_dg, std_error


parser = argparse.ArgumentParser()
parser.add_argument('--df')
parser.add_argument('-d', "--delim", help="Give the csv files delimiter")
parser.add_argument('-n', "--name", help="The id column name")
args = parser.parse_args()


# get working dir
wd = os.getcwd()
print(wd)

df = pd.read_csv(args.df, delimiter=args.delim)
print("HELLO")

mu_exs = []
mu_ex_errs = []
for i in range(len(df)):
    x = df.iloc[i]
    id = str(x[args.name])
    print(id)
    try:
        mu, err = calc_mu_ex(id)
    except FileNotFoundError:
        mu, err = np.nan, np.nan
    mu_exs.append(mu)
    mu_ex_errs.append(err)

df["Restraints"] = mu_exs
df["Restraints_err"] = mu_ex_errs

df.to_csv("restraints_on.txt", index=False)

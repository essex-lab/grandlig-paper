"""
Simple script to take in complex leg, muex, 
restraints leg sym numbers etc and calculate a final free enrgy using an ABFE thermodynamic cycle. 

Where multiple bidning mdoes exisit, a boltzmann average is taken
"""
# Modules to import
import argparse
import numpy as np
import pandas as pd
from openmm.unit import *

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("--comp", help="complex dGs", type=str)
parser.add_argument("--solv", help="solv dGs", type=str)
parser.add_argument("--reston", help="rest on dGs", type=str)
parser.add_argument("--restoff", help="rest off dGs", type=str)
parser.add_argument("--ligands", help="Unique ligands", type=str)

parser.add_argument("--sym", default=None, help="complex dGs", type=str)
args = parser.parse_args()


kT = 298 * kelvin * BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA
kT = kT.in_units_of(kilocalories_per_mole)
def jacobian(dgs, errs, kT=kT):
    dgs = np.asarray(dgs)
    delta_dgs = np.asarray(errs)
    beta = 1 / kT

    exp_factors = np.exp(-beta * dgs)
    sumup = np.sum(exp_factors)  # partition factor
    w = exp_factors / sumup  # weights

    dg_all = (-1 / beta) * np.log(sumup)

    # Error propagation for the sum of Boltzmann factors
    sigma_G = np.sqrt(np.sum((w**2) * (delta_dgs**2)))
    sigma_G_ = np.sqrt(np.sum(((w * delta_dgs) / sumup) ** 2))

    print(sigma_G, sigma_G_)
    return [dg_all._value, sigma_G]


complex = pd.read_csv(args.comp)
solv = pd.read_csv(args.solv)
solv.rename(columns={"Mu_Ex": "Solv", "Mu_Ex_err": "Solv_err"}, inplace=True)

reston = pd.read_csv(args.reston)
restoff = pd.read_csv(args.restoff)
ligands = pd.read_csv(args.ligands)
sym_number = pd.read_csv("sym_nums.txt")

sym_number["sym_corr"] = - kT._value * np.log(sym_number["sym"])

ligands = complex["Ligand"]

unique_ligands = pd.read_csv(
    "/home/will/data_6/LIG_GCNCMC_PAPER/REVIEW_ANALYSIS/fep/unique_ligs.txt"
)
unique_ligands = unique_ligands[["Ligand"]]
unique_ligands.rename(columns={"Ligand": "Name"}, inplace=True)


all_data = complex.merge(solv).merge(reston).merge(restoff).merge(sym_number)


all_data["dG"] = (
    -all_data["Solv"]
    - all_data["Rest_off"]
    + all_data["Complex"]
    - all_data["Restraints"]
)
all_data["dG_err"] = np.sqrt(
    all_data["Complex_err"] ** 2
    + all_data["Solv_err"] ** 2
    + all_data["Restraints_err"] ** 2,
)

all_data["dG_sym"] = all_data["dG"] + all_data["sym_corr"]
all_data["Ligand_stripped"] = all_data["Ligand"].str.split("_").str[0]
jacos_sym = [
    jacobian(
        all_data[all_data["Ligand_stripped"] == i][
            "dG_sym"
        ].tolist(),
        all_data[all_data["Ligand_stripped"] == i][
            "dG_err"
        ].tolist(),
    )
    for i in unique_ligands["Name"]
]

most_fav = [
    np.min(all_data[all_data["Ligand_stripped"] == i]["dG_sym"].tolist())
    for i in unique_ligands["Name"]
]
most_fav_ids = [
    np.argmin(all_data[all_data["Ligand_stripped"] == i]["dG_sym"].tolist())
    for i in unique_ligands["Name"]
]
most_fav_err = [
    all_data[all_data["Ligand_stripped"] == i]["dG_err"].tolist()[most_fav_ids[j]]
    for j, i in enumerate(unique_ligands["Name"])]

jacos_sym = np.asarray(jacos_sym)


unique_ligands["final_dG_sym"] = jacos_sym[:, 0]
unique_ligands["final_err"] = jacos_sym[:, 1]


unique_ligands["most_fav"] = most_fav
unique_ligands["most_fav_err"] = most_fav_err

# unique_ligands.to_csv(
#     "final_FEP.txt",
#     sep="\t",
#     columns=["final_dG", "final_err", "Name"],
#     header=False,
#     index=False,
# )
unique_ligands.to_csv(
    "final_FEP_sym.txt",
    columns=["final_dG_sym", "final_err", "Name"],
    sep="\t",
    header=False,
    index=False,
)
unique_ligands.to_csv(
    "final_FEP_mostfav.txt",
    columns=["most_fav", "most_fav_err", "Name"],
    sep="\t",
    header=False,
    index=False,
)

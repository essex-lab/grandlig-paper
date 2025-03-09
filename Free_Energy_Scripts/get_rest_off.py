import pandas as pd
import argparse
from pymbar.timeseries import statistical_inefficiency, subsample_correlated_data, detect_equilibration
from pymbar import MBAR
import numpy as np
import os
from openmm.unit import *

def calc_mu_ex(id, n_repeats=3):
    WORK_DIR = os.getcwd()
    ligand = id
    #print(WORK_DIR)
    #print(ligand)
    repeats = range(1, n_repeats+1)


    dg_repeats = []
    u_kln_list = []
    for repeat in repeats:
        if not os.path.isfile(f"{ligand}/repeat_{repeat}/lambda_0/dG_log_file.txt"):
            continue
        with open(f"{ligand}/repeat_{repeat}/lambda_0/dG_log_file.txt", 'r') as f:
            for line in f.readlines():
                if line.startswith("The computed standard state correction is"):
                    dg = line.split()[-2]
        return dg


parser = argparse.ArgumentParser()
parser.add_argument('--df')
parser.add_argument('-d', "--delim", help="Give the csv files delimiter")
parser.add_argument('-n', "--name", help="The id column name")
args = parser.parse_args()


# get working dir
wd = os.getcwd()
print(wd)

df = pd.read_csv(args.df, delimiter=args.delim)

mu_exs = []
mu_ex_errs = []
for i in range(len(df)):
    x = df.iloc[i]
    id = str(x[args.name])
    print(id)
    dg = calc_mu_ex(id)
    mu_exs.append(dg)

df["Rest_off"] = mu_exs

df.to_csv("Restraint_off.txt", index=False)


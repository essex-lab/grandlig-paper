# Modules to import
import argparse
from sys import stdout
import numpy as np
from simtk.unit import *
import pymbar
import MDAnalysis as mda
import MDAnalysis.analysis.distances as mda_distances
import MDAnalysis.core.topologyobjects as mda_topo
import mdtraj as md
from matplotlib import pyplot as plt



# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gs", help="dGs", default=5.0, nargs='+', type=float)
parser.add_argument("-T", "--temp", help="Temperature", default=298.0)

args = parser.parse_args()

temp = args.temp * kelvin
# dGs = [g * kilocalories_per_mole for g in args.gs]

kT = BOLTZMANN_CONSTANT_kB * AVOGADRO_CONSTANT_NA * temp
beta = 1/kT
#
# exp_beta_dGs = [np.exp(-beta * g) for g in dGs]
#
# total_dG  = -1/beta * np.log(np.sum(exp_beta_dGs))
# print(total_dG.in_units_of(kilocalories_per_mole))

ligs = {}


with open("primary_fes.txt", 'r') as fi:
    for line in fi.readlines():
        split = line.split()
        ligand = split[0].split("_")[0]
        if ligand not in list(ligs.keys()):
            ligs[ligand] = [float(split[1]) * kilocalories_per_mole]

with open("secondary_fes.txt", 'r') as fi:
    for line in fi.readlines():
        split = line.split()
        ligand = split[0].split("_")[0]
        if ligand not in list(ligs.keys()):
            ligs[ligand] = [float(split[1]) * kilocalories_per_mole]
        else:
            ligs[ligand].append(float(split[1]) * kilocalories_per_mole)


jacobed = []
for l in ligs.keys():
    dGs = ligs[l]
    exp_beta_dGs = [np.exp(-beta * g) for g in dGs]
    total_dG = -1 / beta * np.log(np.sum(exp_beta_dGs))
    jacobed.append(total_dG.value_in_unit(kilocalories_per_mole))

with open("../jacobianed_complex.txt", "w") as fo:
    for i in range(len(ligs.keys())):
        fo.write(f"{jacobed[i]}\t 0.00\t {list(ligs.keys())[i]}\n")


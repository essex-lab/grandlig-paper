from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
import numpy as np
import openmmtools
import pandas as pd
from openff.toolkit import Topology, Molecule

import grandlig as grand

import argparse


parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdb", type=str)
parser.add_argument("-l", "--lig", type=str)
args = parser.parse_args()


pdb = PDBFile(args.pdb)
omm_top = pdb.topology
lig = Molecule.from_file(args.lig)

#top = Topology.from_pdb(args.pdb, unique_molecules=[lig])


#omm_top = top.to_openmm()
# Create system
ff = ForceField("amber14/tip3p.xml")

# Create the SMIRNOFF template generator with the default installed force field (openff-2.1.0)
from openmmforcefields.generators import (
    GAFFTemplateGenerator,
)

gaff = GAFFTemplateGenerator(molecules=lig, forcefield="gaff-1.8")

ff.registerTemplateGenerator(gaff.generator)

list_of_resis = (
    []
)  # Get a list of resids so we can choose a random one to decouple
resname = "L02"

for residue in omm_top.residues():
    if residue.name == resname:
        list_of_resis.append(residue.index)

resid = np.random.choice(list_of_resis)


# Create system

system = ff.createSystem(
    omm_top,
    nonbondedMethod=PME,
    nonbondedCutoff=12.0 * angstroms,
    switchDistance=10.0 * angstroms,
    constraints=HBonds,
)


# Run free energy calculation using grand
grand.potential.calc_mu_ex(system=system,
    topology=omm_top,
    positions=pdb.positions,
    box_vectors=omm_top.getPeriodicBoxVectors(),
    temperature=298.15 * kelvin,
    pressure=1 * bar,
    resname=resname,
    resid=resid,
    n_lambdas=20,
    n_samples=5000,
    n_equil=400,
    log_file="dG.txt",
)
import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
import grandlig as grand
from sys import stdout
import glob

# Load in PDB
# Sphere data
ref_atoms = [
    {"name": "C34", "resname": "MGO", "resid": "1"},
    {"name": "C7", "resname": "MGO", "resid": "1"},
]

sphere_rad = 5.0 * angstrom


# Setup the output files
# Move ghost waters out of the simulation cell
trj = grand.utils.shift_ghost_waters(
    ghost_file="ncmc-ghost-ligs.txt",
    topology="betaCDWithghosts.pdb",
    trajectory="betaCD_raw.dcd",
    output="betaCD_shifted.dcd",
)


# Write out a trajectory of the GCMC sphere
grand.utils.write_sphere_traj(
    radius=sphere_rad,
    ref_atoms=ref_atoms,
    topology="betaCDWithghosts.pdb",
    trajectory="betaCD_raw.dcd",
    output="gcmc_sphere.pdb",
    initial_frame=True,
)

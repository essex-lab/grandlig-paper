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
    {"name": "CG", "resname": "LEU", "resid": "85"},
    {"name": "CB", "resname": "ALA", "resid": "100"},
]

sphere_rad = 8.0 * angstrom


# Setup the output files
# Move ghost waters out of the simulation cell
trj = grand.utils.shift_ghost_molecules(
    ghost_file="ncmc-ghost-ligs.txt",
    topology="T4L99AWithghosts.pdb",
    trajectory="T4L99A_raw.dcd",
)


# Align the trajectory to the protein
grand.utils.align_traj(t=trj, output="T4NCMCholo.dcd", reference="T4L99AWithghosts.pdb")

# Write out a trajectory of the GCMC sphere
grand.utils.write_sphere_traj(
    radius=sphere_rad,
    ref_atoms=ref_atoms,
    topology="T4L99AWithghosts.pdb",
    trajectory="T4NCMCholo.dcd",
    output="gcmc_sphere.pdb",
    initial_frame=True,
)

grand.utils.cluster_molecules(
    topology="T4Withghosts.pdb",
    trajectory="T4NCMCholo.dcd",
    resname="LIG",
    sphere_radius=sphere_rad,
    ref_atoms=ref_atoms,
    cutoff=2.4,
    output="gcmc-clusters.pdb",
)

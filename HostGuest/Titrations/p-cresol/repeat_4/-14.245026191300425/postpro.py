import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
import grandlig as grand
from sys import stdout
import glob

# Move ghost waters out of the simulation cell
trj = grand.utils.shift_ghost_waters(
    ghost_file="ncmc-ghost-ligs.txt",
    topology="betaCDWithghosts.pdb",
    trajectory="betaCD_aligned.dcd",
    output="betaCD_shifted.dcd",
)

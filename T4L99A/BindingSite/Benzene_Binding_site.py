"""
Script that uses GCNCMC to find the occluded Benzene binding site of T4L99A.
"""

import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
import grand
from sys import stdout
import glob

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdb", help="Input PDB file of an equilibrated protein solvent system", type=str)
parser.add_argument("-l", "--lig", help="Input ligand PDB file", type=str)
parser.add_argument("-x", "--xml", help="Input ligand XML file", type=str)
parser.add_argument("-r", "--resname", help="Resname of ligand.", type=str, default='L02')
parser.add_argument("-c", "--concentration", help="Input the concentration you would like to simulate at.", type=float,
                    default=0.5)
args = parser.parse_args()


# Set MD and GCNCMC variables
nonbondedMethod = PME
nonbondedCutoff = 1.2*nanometers
switchDistance=10.0*angstroms
ewaldErrorTolerance = 0.0005  # Default
constraints = HBonds
rigidWater = True  # Default
hydrogenMass = 1*amu

# Integration Options
dt = 0.002*picoseconds
temperature = 298*kelvin
friction = 1.0/picosecond

# GCNCMC params
mu_ex = -0.679*kilocalories_per_mole
concentration = args.conc * molar
volume_per_lig = grand.utils.convert_conc_to_volume(concentration)
n_pert = 1499
n_prop = 50
adams_value = None

# Load in equilibrated protein+solvent PDB
pdb = PDBFile(args.pdb)


# Add ghost molecules
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology,
                                                             pdb.positions,
                                                             molfile=args.lig,
                                                             n=100,
                                                             pdb='T4L99AWithghosts.pdb')

# Load force field and create system
ff = ForceField('amber14-all.xml', 'amber14/tip3p.xml', args.xml)
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=PME,
                         nonbondedCutoff=12.0*angstroms,
                         switchDistance=10.0*angstroms,
                         constraints=HBonds)

# Langevin Integrator
integrator = BAOABIntegrator(temperature, friction, dt)

# Sphere data
ref_atoms = [{'name': 'CA', 'resname': 'PHE', 'resid': '105'},
             {'name': 'CA', 'resname': 'GLU', 'resid': '12'}]
sphere_rad = 26.5 * angstrom


# Set up NCMC Sphere sampler
print('Setting up sphere sampler object...')
ncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(system=system, topology=pdb.topology,
                                                            temperature=temperature, integrator=integrator,
                                                            excessChemicalPotential=mu_ex,
                                                            standardVolume=volume_per_lig,
                                                            referenceAtoms=ref_atoms, sphereRadius=sphere_rad,
                                                            nPertSteps=n_pert, nPropStepsPerPert=n_prop,
                                                            resname=args.resname, ghostFile='ncmc-ghost-ligs.txt',
                                                            log='T4L99A.log', dcd='T4L99A_raw.dcd',
                                                            rst='T4L99A-rst.rst7', overwrite=True,
                                                            recordTraj=False)

# Get platform
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')
print(f'Running on the {platform.getName()} platform.')

# Set up simulation object
simulation = Simulation(pdb.topology, system, ncmc_mover.integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

simulation.reporters.append(StateDataReporter(stdout, 500, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, volume=True))

print('Initialising...')
# Switch off ghosts
ncmc_mover.initialise(simulation.context, simulation, ghosts)
forces = simulation.context.getState(getForces=True).getForces(asNumpy=True)

print('Minimising...')
simulation.minimizeEnergy(tolerance=0.0001*kilojoule/mole, maxIterations=10000)  # Quick Minimisation

print('\nProduction....')
print('Running NCMC...')
for i in range(1000):  #  Number of MD - GCNCMC move cycles to perform
    simulation.step(2000)  # Stop normal MD steps to propagate the system
    ncmc_mover.move(simulation.context, 1)  # Perform a move
    ncmc_mover.report(simulation)  # Report the frame to the DCD file


# Setup the output files
# Move ghost waters out of the simulation cell for better visulisation
trj = grand.utils.shift_ghost_waters(ghost_file='ncmc-ghost-ligs.txt',
                                     topology='T4L99AWithghosts.pdb',
                                     trajectory='T4L99A_raw.dcd')

# Align the trajectory to the protein
grand.utils.align_traj(t=trj, output='T4L99A.dcd', reference="T4L99AWithghosts.pdb")


# Write out a trajectory of the GCMC sphere
grand.utils.write_sphere_traj(radius=sphere_rad,
                              ref_atoms=ref_atoms,
                              topology='T4L99AWithghosts.pdb',
                              trajectory='T4L99A.dcd',
                              output='gcmc_sphere.pdb',
                              initial_frame=True)



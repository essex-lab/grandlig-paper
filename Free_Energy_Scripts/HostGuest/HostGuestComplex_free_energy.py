#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to calculate the free energy accociated with the addition of restraints between a small molecule from its receptor using a
harmonic/boresch restraint to keep the ligand in the pocket
Accepts the system prmtop and inpcrd files, the resname of the small molecule, number of lambda windows and number of samples
- Will Poole
"""

# Modules to import
import argparse
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from sys import stdout
import numpy as np
import openmmtools
from openmmtools.constants import STANDARD_STATE_VOLUME
import pymbar
import MDAnalysis as mda
import MDAnalysis.analysis.distances as mda_distances
import MDAnalysis.core.topologyobjects as mda_topo
import mdtraj as md
import grandlig as grand

# Functions
def get_lambda_values(lambda_in, ratio='half'):
    """
    Calculate the lambda_sterics and lambda_electrostatics values for a given lambda.
    Electrostatics are decoupled from lambda=1 to 0.5, and sterics are decoupled from
    lambda=0.5 to 0.
    Parameters
    ----------
    lambda_in : float
        Input lambda value
    Returns
    -------
    lambda_vdw : float
        Lambda value for steric interactions
    lambda_ele : float
        Lambda value for electrostatic interactions
    """
    if ratio == 'half':
        x = 2.0
        y = 2.0
        z = 0.5
    if ratio == 'quarter':
        x = 1.325
        y = 4.00
        z = 0.75
    if lambda_in > 1.0:
        # Set both values to 1.0 if lambda > 1
        lambda_vdw = 1.0
        lambda_ele = 1.0
    elif lambda_in < 0.0:
        # Set both values to 0.0 if lambda < 0
        lambda_vdw = 0.0
        lambda_ele = 0.0
    else:
        # Scale values between 0 and 1
        lambda_vdw = min([1.0, x*lambda_in])
        lambda_ele = max([0.0, y*(lambda_in-z)])
    return lambda_vdw, lambda_ele


def change_lambda(context, new_lambda, rest):
    rest_lam = rest[new_lambda]
    context.setParameter('lambda_restraints', rest_lam)
    return None

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdb", help="Input PDB file of the water + small molecule system", default='system1.pdb')
parser.add_argument("-x", "--xml", help="Input ligand XML file File", default='system.dcd')
parser.add_argument("-r", "--resname", help="Resname of ligand.", type=str, default='LIG')
parser.add_argument("-lam", "--lam", help="Input the actual lambda value were at between 0-39", default=0, type=int)
parser.add_argument("-n_sam", "--n_samples", help="Input the number of samples to collect per lambda window.", type=int, default=1000)
parser.add_argument("-T", '--temp', help='Input simulation temperature in kelvin.', type=float, default=298.0)
args = parser.parse_args()

# Set up variables
resname = args.resname

n_samples = args.n_samples
lam_index = args.lam
n_steps = 1500
prefix = args.pdb.split('.pdb')[0][-1]
log_file = 'dG_log_file.txt'

# Set up lambdas and Energy matrix
n_lambdas = 40
eleclambdas = [x for x in np.linspace(1.0, 0.0, 10)] + [0] * (n_lambdas - 10)
vdwlambdas = [1] * 9 + [x for x in np.linspace(1.0, 0.1, (n_lambdas - 19))] + [0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.00]

lambdas = np.linspace(1.0, 0.0, n_lambdas)
U = np.zeros((n_lambdas, n_lambdas, n_samples))

ele = eleclambdas[lam_index]
vdw = vdwlambdas[lam_index]

with open(log_file, 'w') as log:
    log.write(f'This the log file for the free energy calculation of resdiue: {resname} in complex with ERK2.\n')
    log.write(f'Simulating {n_lambdas} lambda windows and collecting {n_samples} samples per lambda. '
              f'Simulations are being run for {n_steps} steps ({n_steps * 0.002} ps)\n')

# Set up basic OpenMM system
pdb = PDBFile(args.pdb)

# Setup restraint parameters
ref_atoms = [{'name': 'C34', 'resname': 'MGO', 'resid': '1'},
             {'name': 'C7', 'resname': 'MGO', 'resid': '1'}]

protein_group = []
for point in ref_atoms:
    for atom in pdb.topology.atoms():
        if atom.residue.id == point['resid']:
            if atom.name == point['name']:
                protein_group.append(int(atom.index))

print(protein_group)

ligand_group = []
for residue in pdb.topology.residues():
    if residue.name == resname:
        for atom in residue.atoms():
            print(atom)
            ligand_group.append(int(atom.index))

print(ligand_group)

for residue in pdb.topology.residues():
    #print(residue)
    if residue.name == resname:
        resid = residue.index
     #   print(residue, resid, residue.id)


pressure = 1*bar
temperature = args.temp*kelvin
integrator = openmmtools.integrators.BAOABIntegrator(temperature, 1/picosecond, 0.002*picoseconds)
barostat = MonteCarloBarostat(pressure, temperature, 25)
ff = ForceField('../../../bCD_FINAL.xml', 'amber14/tip3p.xml', args.xml)
topology = pdb.topology
system = ff.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.2*nanometer, constraints=HBonds,
                             switchDistance=10*angstroms)

kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * integrator.getTemperature()  # Define kT

# Define a GCMC sampler object, just to allow easy switching of a water - won't use this to sample
gcmc_mover = grand.samplers.BaseGrandCanonicalMonteCarloSampler(system=system, topology=topology,
                                                                temperature=temperature, resname=resname,
                                                                log=log_file,
                                                                ghostFile='calc_mu-ghosts.txt',
                                                                overwrite=True)

restraint_force = openmmtools.forces.FlatBottomRestraintForce(0.6 * (kilocalories_per_mole / (angstroms)**2),
                                                              5 * angstroms,
                                                              protein_group, ligand_group, 'lambda_restraints')
system.addForce(restraint_force)

if pressure is not None:
    system.addForce(MonteCarloBarostat(pressure, temperature, 25))

thermo_state = openmmtools.states.ThermodynamicState(system=system, temperature=temperature, pressure=pressure)

# Define platform and set precision - try CUDA first, then OpenCL - (CUDA will be fine for me)
try:
    platform = Platform.getPlatformByName('CUDA')
    platform.setPropertyDefaultValue('Precision', 'mixed')
    print('Running on CUDA platform...')
except:
    try:
        platform = Platform.getPlatformByName('OpenCL')
        platform.setPropertyDefaultValue('Precision', 'mixed')
        print('Running on OpenCL platform...')
    except:
        raise Exception('No CUDA or OpenCL platform found!')

# Initialise Simulation object
simulation = Simulation(topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature)
original_box_vectors = pdb.topology.getPeriodicBoxVectors()
simulation.context.setPeriodicBoxVectors(*original_box_vectors)

print('Simulation Created')
# Make sure the GCMC sampler has access to the Context
gcmc_mover.context = simulation.context

initial_positions = simulation.context.getState(getPositions=True).getPositions()  # Get the initial positions
with open(f"{prefix}_openMM.pdb", 'w') as f:  # Write out the initial structure in OpenMM coords - useful
    pdb.writeFile(positions=initial_positions, topology=simulation.topology, file=f)

#simulation.reporters.append(DCDReporter('simulation.dcd', 10000))  # Write a frame every 20 ps for checking purposes
simulation.reporters.append(StateDataReporter('simulation_log.txt', 1000, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, volume=True))

#with open(log_file, 'a') as log:
#    log.write(f'The computed standard state correction is {(std_state).in_units_of(kilocalories_per_mole)}\n')

#simulation.minimizeEnergy(tolerance=0.0001*kilojoule/mole, maxIterations=10000)  # Quick Minimisation

i = lam_index

# Commence simulation
with open(log_file, 'a') as log:
    # Set lambda
    log.write('Setting Restraint lambda..............\n')
    change_lambda(simulation.context, 0, [1])
    log.write('Simulating with lambda_restraint = {}\n'.format(simulation.context.getParameter('lambda_restraints')))

# Simulate the system at specified lambda window
# Set lambda values
print('Simulating at lambda index = {:.4f}, elec = {}, vdw = {}'.format(np.round(lambdas[i], 4), ele, vdw))
gcmc_mover.logger.info('Simulating at lambda = {:.4f}'.format(np.round(lambdas[i], 4)))
gcmc_mover.adjustSpecificMolecule(resid, lambdas[i], ele=ele, vdw=vdw)

# Minimise
simulation.minimizeEnergy()
# Equilibrate

n_steps = (3 * nanosecond) / (0.002 * picosecond)
print('Equilibrating for {} steps at lambda = {}'.format(n_steps, np.round(lambdas[i], 4)))

simulation.step(int(n_steps))
print('Equil Done.. Simulation now')
n_equil = 1500

for k in range(n_samples):
    # Run production MD
    simulation.step(n_equil)
    box_vectors = simulation.context.getState(getPositions=True).getPeriodicBoxVectors()
    volume = box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]
    # Calculate energy at each lambda value
    for j in range(n_lambdas):
        # Set lambda value
        gcmc_mover.adjustSpecificMolecule(resid, lambdas[j], ele=eleclambdas[j], vdw=vdwlambdas[j])
        # Calculate energy
        U[i, j, k] = simulation.context.getState(getEnergy=True).getPotentialEnergy() / gcmc_mover.kT
    # Reset lambda value
    gcmc_mover.adjustSpecificMolecule(resid, lambdas[i], ele=eleclambdas[i], vdw=vdwlambdas[i])

# Save the numpy matrix (for now)
np.save(f'U_matrix_{i}.npy', U)

# Calculate equilibration & number of uncorrelated samples
N_k = np.zeros(n_lambdas, np.int32)
for i in range(n_lambdas):
    n_equil, g, neff_max = pymbar.timeseries.detectEquilibration(U[i, i, :])
    indices = pymbar.timeseries.subsampleCorrelatedData(U[i, i, :], g=g)
    N_k[i] = len(indices)
    U[i, :, 0:N_k[i]] = U[i, :, indices].T


gcmc_mover.logger.info(f'Std State = {restraint_force.compute_standard_state_correction(thermodynamic_state=thermo_state, square_well=True, radius_cutoff=5*angstroms)}\n')

# Calculate free energy differences
mbar = pymbar.MBAR(U, N_k)
results = mbar.getFreeEnergyDifferences()
deltaG_ij = results[0]
ddeltaG_ij = results[1]

# Extract overall free energy change
dG = deltaG_ij[0, -1]

# Write out intermediate free energies
for i in range(n_lambdas):
    dG_i = (deltaG_ij[0, i] * gcmc_mover.kT).in_units_of(kilocalorie_per_mole)
    gcmc_mover.logger.info('Free energy ({:.3f} -> {:.3f}) = {}'.format(lambdas[0], lambdas[i], dG_i))

# Convert free energy to kcal/mol
dG = (dG * gcmc_mover.kT).in_units_of(kilocalorie_per_mole)
dG_err = (ddeltaG_ij[0, -1] * gcmc_mover.kT).in_units_of(kilocalorie_per_mole)

gcmc_mover.logger.info('Excess chemical potential = {}'.format(dG))
gcmc_mover.logger.info('Estimated error = {}'.format(dG_err))
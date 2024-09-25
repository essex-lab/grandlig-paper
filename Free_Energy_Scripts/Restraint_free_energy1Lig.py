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
import numpy as np
from openmmtools.integrators import BAOABIntegrator
import pymbar
from boresch_rest import BoreschRestraint
import grandlig as grand

def change_lambda(context, new_lambda, rest):
    rest_lam = rest[new_lambda]
    context.setParameter('lambda_boresch', rest_lam)
    return None

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdb", help="Input PDB file of the water + small molecule system", default='system1.pdb')
parser.add_argument("-x", "--xml", help="Input ligand XML file File", default='system.dcd')
parser.add_argument("-r", "--resname", help="Resname of ligand.", type=str, default='LIG')
parser.add_argument("-a", "--restrained_atoms", help="Atoms to be restrained with a Boresch restraint - ZERO INDEXED", nargs='+', type=int, default=None)
parser.add_argument("-ref", "--ref_values", help="Reference values for the restraints r a1 a2 d1 d2 d3", type=float, nargs="+", default=None)
parser.add_argument("-lam", "--lam", help="Input the actual lambda value were at between 0-29", default=0, type=int)
parser.add_argument("-n_sam", "--n_samples", help="Input the number of samples to collect per lambda window.", type=int, default=1000)
parser.add_argument("-T", '--temp', help='Input simulation temperature in kelvin.', type=float, default=298.0)
parser.add_argument("-K", '--rest_str', help='Input strength of the restraint i kcal/mol.', type=float, default=10)

args = parser.parse_args()

# Set up variables
resname = args.resname
n_samples = args.n_samples
lam_index = args.lam
n_steps = 1500
prefix = args.pdb.split('.pdb')[0][-1]
log_file = 'dG_log_file.txt'

# Set up lambdas and Energy matrix
restraint_lambdas = [0.00, 0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.85, 1.0]
n_lambdas = len(restraint_lambdas)
print(n_lambdas)
lambdas = np.linspace(0.0, 1.0, n_lambdas)
U = np.zeros((n_lambdas, n_lambdas, n_samples))

with open(log_file, 'w') as log:
    log.write(f'This the log file for the free energy calculation of resdiue: {resname} in complex with protein.\n')
    log.write(f'Simulating {n_lambdas} lambda windows and collecting {n_samples} samples per lambda. This simulation is for lambda = {restraint_lambdas[lam_index]} '
              f'Simulations are being run for {n_steps} steps ({n_steps * 0.002} ps)\n')

# Set up basic OpenMM system
pdb = PDBFile(args.pdb)
pressure = 1*bar
temperature = args.temp*kelvin
integrator = BAOABIntegrator(temperature, 1/picosecond, 0.002*picoseconds)
barostat = MonteCarloBarostat(pressure, temperature, 25)
ff = ForceField('amber14-all.xml', 'amber14/tip3p.xml', args.xml)
topology = pdb.topology
system = ff.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1.2*nanometer, constraints=HBonds,
                             switchDistance=10*angstroms)


kT = AVOGADRO_CONSTANT_NA * BOLTZMANN_CONSTANT_kB * integrator.getTemperature()  # Define kT

# Define a GCMC sampler object, just to allow easy switching of a water - won't use this to sample
gcmc_mover = grand.samplers.BaseGrandCanonicalMonteCarloSampler(system=system, topology=topology,
                                                                temperature=temperature, resname=args.resname,
                                                                log=log_file,
                                                                ghostFile='calc_mu-ghosts.txt',
                                                                overwrite=True)

restrained_atoms = args.restrained_atoms
print(restrained_atoms)

if args.ref_values == None:
    r_ref = None
    a_ref = None
    d_ref = None
else:
    r_ref = args.ref_values[0]
    a_ref = args.ref_values[1:3]
    d_ref = args.ref_values[3:]


rest_str = args.rest_str

boresch_restraints, std_state = BoreschRestraint(pdb, restrained_atoms=restrained_atoms,
                                                 K_r=rest_str * kilocalories_per_mole / angstroms**2,
                                                 K_theta=rest_str * kilocalories_per_mole / radians**2,
                                                 K_phi=rest_str * kilocalories_per_mole / radians**2,
                                                 r=r_ref, a=a_ref, d=d_ref, lambda_controlled='lambda_boresch', kT=kT)
system.addForce(boresch_restraints)
if pressure is not None:
    system.addForce(MonteCarloBarostat(pressure, temperature, 25))

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

simulation.reporters.append(DCDReporter('simulation.dcd', 10000))  # Write a frame every 20 ps for checking purposes
simulation.reporters.append(StateDataReporter('simulation_log.txt', 1000, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, volume=True))

with open(log_file, 'a') as log:
    log.write(f'The computed standard state correction is {(std_state).in_units_of(kilocalories_per_mole)}\n')

i = lam_index

change_lambda(simulation.context, 0, restraint_lambdas)
simulation.minimizeEnergy(tolerance=0.00001*kilojoule/mole, maxIterations=10000)  # Quick Minimisation

change_lambda(simulation.context, i, restraint_lambdas)
simulation.minimizeEnergy(tolerance=0.00001*kilojoule/mole, maxIterations=10000)  # Quick Minimisation

# Commence simulation
with open(log_file, 'a') as log:
    log.write(f'Simulating at lambda value {i}\n')
# Set lambda
    log.write('Setting lambda..............\n')
change_lambda(simulation.context, i, restraint_lambdas)

with open(log_file, 'a') as log:
    log.write('Simulating with lambda_boresch = {}\n'.format(simulation.context.getParameter('lambda_boresch')))

# Minimise
print("Minimize!!")
simulation.minimizeEnergy(tolerance=0.00001*kilojoule/mole, maxIterations=1000000)
# Equilibrate

scale = np.linspace(0, restraint_lambdas[i], 100)
for _ in range(100):
    simulation.context.setParameter('lambda_boresch', scale[_])
    simulation.step(2500)
change_lambda(simulation.context, i, restraint_lambdas)



n_steps = (2 * nanosecond) / (0.002 * picosecond)
print('Equilibrating for {} steps at lambda = {}'.format(n_steps, np.round(restraint_lambdas[i], 4)))

simulation.step(int(n_steps))
print('Equil Done.. Simulation now')

# 1500 * 1000 * 0.002 = 3ns

for k in range(n_samples):
    simulation.step(1500)
    box_vectors = simulation.context.getState(getPositions=True).getPeriodicBoxVectors()
    volume = box_vectors[0][0] * box_vectors[1][1] * box_vectors[2][2]
    # Calculate energy at each lambda value
    for j in range(n_lambdas):
        change_lambda(simulation.context, j, restraint_lambdas)

        # Calculate Energy at eahch lambda state with volume correction?
        U[i, j, k] = (simulation.context.getState(getEnergy=True).getPotentialEnergy() + (pressure * volume * AVOGADRO_CONSTANT_NA) ) / kT
    # Reset lambda value
#    change_lambda(simulation.context, harmonic_restraints, i, restraint_lambdas)
    change_lambda(simulation.context, i, restraint_lambdas)

final_positions = simulation.context.getState(getPositions=True).getPositions()  # Get the initial positions
with open(f"./FinalFrame.pdb", 'w') as f:  # Write out the initial structure in OpenMM coords - useful
    pdb.writeFile(positions=final_positions, topology=simulation.topology, file=f)
# Save Energy numpy array
np.save('U_matrix.npy', U)


# with open(log_file, 'a') as log:
#     log.write('Simulation complete. Performing analysis...\n')
# # Perform Pymbar analysis - copied from Marley need to look
# N_k = np.zeros(n_lambdas, np.int32)
# for i in range(n_lambdas):
#     n_equil, g, neff_max = pymbar.timeseries.detectEquilibration(U[i, i, :])  # Detect equilibrated states for this lambda
#     indicies = pymbar.timeseries.subsampleCorrelatedData(U[i, i, :], g=g)  # identify uncorrelated data
#     N_k[i] = len(indicies)
#     U[i, :, 0:N_k[i]] = U[i, :, indicies].T
#
# # Calculate free energy differences
# mbar = pymbar.MBAR(U, N_k)
# [deltaG_ij, ddeltaG_ij, theta_ij] = mbar.getFreeEnergyDifferences(return_theta=True)
# dG = deltaG_ij[0, -1]
#
# # Write out the intermediate free energies
# with open(log_file, 'a') as log:
#     log.write('Intermediate free energies')
# for i in range(n_lambdas):
#     dG_i = (deltaG_ij[0, i] * kT).in_units_of(kilocalorie_per_mole)
#     with open(log_file, 'a') as log:
#         log.write('Free energy ({:.3f} -> {:.3f}) = {}\n'.format(lambdas[0], lambdas[i], dG_i))
#
# # Convert free energy and error to kcal / mol
# dG = (dG * kT).in_units_of(kilocalorie_per_mole)
# dG_err = (ddeltaG_ij[0, -1] * kT).in_units_of(kilocalorie_per_mole)
#
# with open(log_file, 'a') as log:
#     log.write(f'Restraint free energy = {dG}\n')
#     log.write(f'Estimated error = {dG_err}\n')



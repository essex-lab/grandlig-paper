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
import pymbar
from boresch_rest import BoreschRestraint
import grandlig as grand
from openmmtools.integrators import BAOABIntegrator


# Functions

def change_lambda(context, new_lambda, rest):
    rest_lam = rest[new_lambda]
    context.setParameter('lambda_boresch', rest_lam)
    return None

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--pdb", help="Input PDB file of the water + small molecule system", default='system1.pdb')
parser.add_argument("-x", "--xml", help="Input ligand XML file File", default='system.dcd')
parser.add_argument("-r", "--resname", help="Resname of ligand.", type=str, default='L02')
parser.add_argument("-a", "--restrained_atoms", help="Atoms to be restrained with a Boresch restraint ", nargs='+', type=int, default=None)
parser.add_argument("-ref", "--ref_values", help="Reference values for the restraints r a1 a2 d1 d2 d3", type=float, nargs="+", default=None)
parser.add_argument("-lam", "--lam", help="Input the actual lambda value were at between 0-39", default=0, type=int)
parser.add_argument("-n_sam", "--n_samples", help="Input the number of samples to collect per lambda window.", type=int, default=2000)
parser.add_argument("-T", '--temp', help='Input simulation temperature in kelvin.', type=float, default=298.0)
parser.add_argument("-K", '--rest_str', help='Input strength of the restraint i kcal/mol.', type=float, default=10)
args = parser.parse_args()

# Set up variables
resname = args.resname[0]
n_samples = args.n_samples
lam_index = args.lam
n_equil = 1250  # 1250 * 2000 * 0.002 = 5 ns
prefix = args.pdb.split('.pdb')[0][-1]
log_file = 'dG_log_file.txt'

# Set up lambdas and Energy matrix
n_lambdas = 40
eleclambdas = [x for x in np.linspace(1.0, 0.0, 11)] + [0] * (n_lambdas - 11)
vdwlambdas = [1] * 10 + [x for x in np.linspace(1.0, 0.1, (n_lambdas - 20))] + [0.09, 0.08, 0.07, 0.06, 0.05, 0.04, 0.03, 0.02, 0.01, 0.00]


lambdas = np.linspace(1.0, 0.0, n_lambdas)
U = np.zeros((n_lambdas, n_lambdas, n_samples))

ele = np.round(eleclambdas[lam_index], 3)
vdw = np.round(vdwlambdas[lam_index], 3)

with open(log_file, 'w') as log:
    log.write(f'This the log file for the free energy calculation of resdiue: {resname} in complex.\n')
    log.write(f'Simulating {n_lambdas} lambda windows and collecting {n_samples} samples per lambda. '
              f'Samples are recorded after every {n_equil} steps ({n_equil * 0.002} ps)\n')

# Set up basic OpenMM system
pdb = PDBFile(args.pdb)

for residue in pdb.topology.residues():
    if residue.name == args.resname:
        resid = residue.index
        print(residue, resid, residue.id)

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

print(r_ref, a_ref, d_ref)

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

# Minimize before adding restraint
simulation.context.setParameter('lambda_boresch', 0)
simulation.minimizeEnergy(tolerance=0.00001*kilojoule/mole, maxIterations=1000000)


# simulation.reporters.append(DCDReporter('simulation.dcd', 25000))  # Write a frame every 50 ps for checking purposes
simulation.reporters.append(StateDataReporter('simulation_log.txt', 1000, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, volume=True))

with open(log_file, 'a') as log:
    log.write(f'The computed standard state correction is {(std_state).in_units_of(kilocalories_per_mole)}\n')

i = lam_index

# Commence simulation
with open(log_file, 'a') as log:
    # Set lambda
    log.write('Setting Boresch lambda..............\n')
    simulation.context.setParameter('lambda_boresch', 1)
    log.write('Simulating with lambda_boresch = {}\n'.format(simulation.context.getParameter('lambda_boresch')))

# Simulate the system at specified lambda window
# Set lambda values
print('Simulating at lambda index = {:.4f}, elec = {}, vdw = {}'.format(np.round(lambdas[i], 4), ele, vdw))
gcmc_mover.logger.info('Simulating at lambda = {:.4f}'.format(np.round(lambdas[i], 4)))
gcmc_mover.adjustSpecificMolecule(resid, lambdas[i], ele=ele, vdw=vdw)

# Minimise
simulation.minimizeEnergy(tolerance=0.00001*kilojoule/mole, maxIterations=1000000)

# Equilibrate
n_steps = (2 * nanosecond) / (0.002 * picosecond)
print('Equilibrating for {} steps at lambda = {}'.format(n_steps, np.round(lambdas[i], 4)))

simulation.step(int(n_steps))
print('Equil Done.. Simulation now')
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

final_positions = simulation.context.getState(getPositions=True).getPositions()  # Get the initial positions
with open(f"./FinalFrame.pdb", 'w') as f:  # Write out the initial structure in OpenMM coords - useful
    pdb.writeFile(positions=final_positions, topology=simulation.topology, file=f)

# Save the numpy matrix (for now)
np.save(f'U_matrix_{i}.npy', U)

# # Calculate equilibration & number of uncorrelated samples
# N_k = np.zeros(n_lambdas, np.int32)
# for i in range(n_lambdas):
#     n_equil, g, neff_max = pymbar.timeseries.detectEquilibration(U[i, i, :])
#     indices = pymbar.timeseries.subsampleCorrelatedData(U[i, i, :], g=g)
#     N_k[i] = len(indices)
#     U[i, :, 0:N_k[i]] = U[i, :, indices].T
#
# # Calculate free energy differences
# mbar = pymbar.MBAR(U, N_k)
# results = mbar.getFreeEnergyDifferences()
# deltaG_ij = results[0]
# ddeltaG_ij = results[1]
#
# # Extract overall free energy change
# dG = deltaG_ij[0, -1]
#
# # Write out intermediate free energies
# for i in range(n_lambdas):
#     dG_i = (deltaG_ij[0, i] * gcmc_mover.kT).in_units_of(kilocalorie_per_mole)
#     gcmc_mover.logger.info('Free energy ({:.3f} -> {:.3f}) = {}'.format(lambdas[0], lambdas[i], dG_i))
#
# # Convert free energy to kcal/mol
# dG = (dG * gcmc_mover.kT).in_units_of(kilocalorie_per_mole)
# dG_err = (ddeltaG_ij[0, -1] * gcmc_mover.kT).in_units_of(kilocalorie_per_mole)
#
# gcmc_mover.logger.info('Excess chemical potential = {}'.format(dG))
# gcmc_mover.logger.info('Estimated error = {}'.format(dG_err))
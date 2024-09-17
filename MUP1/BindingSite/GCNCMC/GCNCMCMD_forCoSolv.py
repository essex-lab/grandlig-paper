import argparse
import grandlig as grand
from openmm import *
from openmm.app import *
from openmm.unit import *
import openmmtools
import numpy as np
from openmmtools.integrators import BAOABIntegrator
import math

parser = argparse.ArgumentParser()
parser.add_argument('--pdb')
parser.add_argument('--lig')
parser.add_argument('-x', '--xml')
parser.add_argument('--mu', help='Ligand mu_ex in units of kcal/mol', type=float)
parser.add_argument('--conc', help='Concentration in units of M', type=float)
parser.add_argument('--st', help='Switching time in ps', type=int)
parser.add_argument('--sph', help='Give the anchor atoms for the GCNMC Sphere in form name, resname, resid', nargs='+')
parser.add_argument('--rad', help='GCMC Sphere radius in Angstroms', type=float)

parser.add_argument('--hmr', help='Use HMR?', action='store_true')
args = parser.parse_args()

print(args.hmr)
# Basic System Configuration
nonbondedMethod = PME
nonbondedCutoff = 1.2*nanometers
switchDistance = 10.0*angstroms
ewaldErrorTolerance = 0.0005  # Default
constraints = HBonds
rigidWater = True  # Default
hydrogenMass = 1.0*amu

# Integration Options
dt = 0.002*picoseconds
temperature = 298*kelvin
friction = 1.0/picosecond

# HMR?
if args.hmr:
    dt = 0.004 * picoseconds
    hydrogenMass = 2.0 * amu

# GCNCMC params
mu_ex = args.mu * kilocalories_per_mole
conc = args.conc * molar
v_per_lig = grand.utils.convert_conc_to_volume(conc)

switching_time = args.st * picoseconds
n_prop = 50
n_pert = int((switching_time / (n_prop * dt)) - 1)
print(n_pert)


# GCMC Sphere
sph_name, sph_resname, sph_resi = args.sph
ref_atoms = [{'name': sph_name, 'resname': sph_resname, 'resid': sph_resi}]

sphere_rad = args.rad * angstroms

# Input Files
protein_file = args.pdb
lig_file = args.lig
xml_file = args.xml
pdb = PDBFile(protein_file)

# Get platform
platform = Platform.getPlatformByName('CUDA')
platformProperties = {'Precision': 'mixed'}
platform.setPropertyDefaultValue('Precision', 'mixed')
print(f'Running on {platform.getName()} platform.')

# Prepare Simulation and add ghost molecules for GCMC.
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(pdb.topology,
                                                             pdb.positions,
                                                             molfile=lig_file,
                                                             n=200,
                                                             pdb='Protein_Ghosts.pdb')

# Use the ghost resids to get the ligand resname
for residue in pdb.topology.residues():
    if residue.index == ghosts[0]:
        resname = residue.name


# Load force field and create system
ff = ForceField('amber14-all.xml', 'amber14/tip3p.xml', xml_file)
system = ff.createSystem(pdb.topology,
                         nonbondedMethod=nonbondedMethod,
                         nonbondedCutoff=nonbondedCutoff,
                         switchDistance=switchDistance,
                         constraints=constraints,
                         rigidWater=rigidWater,
                         hydrogenMass=hydrogenMass)

# Langevin Integrator
integrator = BAOABIntegrator(temperature, friction, dt)


# Set up NCMC Sphere sampler
print('Setting up sphere sampler object...')
ncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(system=system, topology=pdb.topology,
                                                            temperature=temperature,
                                                            integrator=integrator,
                                                            excessChemicalPotential=mu_ex,
                                                            adams=None, standardVolume=v_per_lig,
                                                            referenceAtoms=ref_atoms, sphereRadius=sphere_rad,
                                                            nPertSteps=n_pert, nPropStepsPerPert=n_prop, timeStep=dt,
                                                            resname=resname,
                                                            ghostFile='ncmc-ghost-ligs.txt', log='gcncmc.log',
                                                            rst="checkpoint.rst7",
                                                            dcd='gcncmc_raw.dcd', overwrite=True, recordTraj=False, maxN=999)

# Set up simulation object
simulation = Simulation(pdb.topology, system, ncmc_mover.integrator, platform, platformProperties=platformProperties)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())


simulation.reporters.append(StateDataReporter('simulation_log.txt', 1000, step=True, time=True,
        potentialEnergy=True, temperature=True, volume=True, elapsedTime=True))

# Switch off ghosts
ncmc_mover.initialise(simulation.context, simulation, ghosts)

print('Minimising...')
simulation.minimizeEnergy(tolerance=0.0001*kilojoule/mole, maxIterations=100000)  # Quick Minimisation


print('Equilibration...')
# Run some MD for equil
equil_time = 2 * nanosecond
equil_md_steps = int((equil_time) / dt)
print(f"Simulating for {equil_md_steps} MD equil steps at {dt} ps ({equil_time})")
for i in range(200):
    simulation.step(equil_md_steps // 200)
    ncmc_mover.move(simulation.context, 1)

ncmc_mover.reset()  # Reset any tracked variables

print("Running GCNCMC/MD")
production_MD_time = 50 * nanosecond
n_moves_wanted = 1000
md_per_move = production_MD_time / n_moves_wanted
steps_per_move = math.floor((md_per_move) / dt)

frame_write_interval = 10 * picoseconds
steps_per_write = math.ceil(steps_per_move / (md_per_move/frame_write_interval))

print('\nProduction....')

print('Running NCMC...')
for i in range(n_moves_wanted):
    for j in range(steps_per_move):
        simulation.step(1)
        if j % steps_per_write == 0:
            ncmc_mover.report(simulation)
    ncmc_mover.move(simulation.context, 1)

ncmc_mover.report(simulation)  # Final report for move 500

# Setup the output files
# Move ghost waters out of the simulation cell
trj = grand.utils.shift_ghost_waters(ghost_file='ncmc-ghost-ligs.txt',
                                     topology='Protein_Ghosts.pdb',
                                     trajectory='gcncmc_raw.dcd')


# Align the trajectory to the protein
grand.utils.align_traj(t=trj, output='gcncmc.dcd', reference="Protein_Ghosts.pdb")

# Write out a trajectory of the GCMC sphere
grand.utils.write_sphere_traj(radius=sphere_rad,
                              ref_atoms=ref_atoms,
                              topology='Protein_Ghosts.pdb',
                              trajectory='gcncmc.dcd',
                              output='gcmc_sphere.pdb',
                              initial_frame=True)






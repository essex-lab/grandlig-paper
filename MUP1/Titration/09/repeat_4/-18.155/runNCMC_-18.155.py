import math

import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
import grandlig as grand
from sys import stdout
import glob

# Input Files
pdb = PDBFile("../../../ProteinEquiled.pdb")  # Could be an arg_parser
lig_pdb = glob.glob("../../*.pdb")
lig_xml = glob.glob("../../*.xml")

# System Configuration
nonbondedMethod = PME
nonbondedCutoff = 1.2 * nanometers
switchDistance = 10.0 * angstroms
ewaldErrorTolerance = 0.0005  # Default
constraints = HBonds
rigidWater = True  # Default
hydrogenMass = 2.0 * amu

# Integration Options
dt = 0.004 * picoseconds
temperature = 298 * kelvin
friction = 1.0 / picosecond
"""
Dont need these for GCMC
pressure = 1.0*atmospheres
barostatInterval = 25
"""

# GCNCMC params
mu_ex = None
volume_per_lig = None
n_pert = 749
n_prop = 50
adams_value = -18.155


# Get platform
platform = Platform.getPlatformByName("CUDA")
platformProperties = {"Precision": "mixed"}
platform.setPropertyDefaultValue("Precision", "mixed")
print(f"Running on {platform} platform.")

# Prepare Simulation and add ghost molecules for GCMC.
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(
    pdb.topology, pdb.positions, molfile=lig_pdb[0], n=20, pdb="MUP1WithGhosts.pdb"
)
# Load force field and create system
ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", lig_xml[0])
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=nonbondedMethod,
    nonbondedCutoff=nonbondedCutoff,
    switchDistance=switchDistance,
    constraints=constraints,
    rigidWater=rigidWater,
    hydrogenMass=hydrogenMass,
)

# Langevin Integrator
integrator = BAOABIntegrator(temperature, friction, 0.002 * picoseconds)

# Sphere data and gcmc stuff
ref_atoms = [
    {"name": "CA", "resname": "PHE", "resid": "57"},
    {"name": "CA", "resname": "LEU", "resid": "106"},
]
sphere_rad = 5.5 * angstrom

# Set up NCMC Sphere sampler
print("Setting up sphere sampler object...")
ncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(
    system=system,
    topology=pdb.topology,
    temperature=temperature,
    integrator=integrator,
    excessChemicalPotential=mu_ex,
    adams=adams_value,
    standardVolume=volume_per_lig,
    referenceAtoms=ref_atoms,
    sphereRadius=sphere_rad,
    nPertSteps=n_pert,
    nPropStepsPerPert=n_prop,
    timeStep=dt,
    resname="L02",
    ghostFile="ncmc-ghost-ligs.txt",
    log="MUP1NCMC.log",
    rst="checkpoint.rst7",
    dcd="MUP1NCMCholo_raw.dcd",
    overwrite=True,
    recordTraj=False,
    maxN=1,
)

# Set up simulation object
simulation = Simulation(
    pdb.topology,
    system,
    ncmc_mover.integrator,
    platform,
    platformProperties=platformProperties,
)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(298 * kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

simulation.reporters.append(
    StateDataReporter(
        stdout,
        1000,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        volume=True,
    )
)

# Switch off ghosts
ncmc_mover.initialise(simulation.context, simulation, ghosts)
ncmc_mover.deleteMoleculesInGCMCSphere()

print("Minimising...")
simulation.minimizeEnergy(
    tolerance=0.0001 * kilojoules_per_mole / nanometer, maxIterations=10000
)  # Quick Minimisation

print("Equilibration...")
# Run some MD for equil
equil_md_steps = math.ceil((1.0 * nanosecond) / (0.002 * picoseconds))
print(f"Simulating for {equil_md_steps} MD equil steps at 2 fs")
simulation.step(equil_md_steps)

integrator.setStepSize(dt)
print(integrator.getStepSize())

moves_since_accepted = 0
n_acc = 0
for i in range(150):
    simulation.step(1000)
    ncmc_mover.move(simulation.context, 1)
    if ncmc_mover.n_accepted == n_acc + 1:
        n_acc += 1
        moves_since_accepted = 0
    elif ncmc_mover.n_accepted == n_acc:
        moves_since_accepted += 1


print(
    "{}/{} equilibration NCMC moves accepted. N = {}".format(
        ncmc_mover.tracked_variables["n_accepted"],
        ncmc_mover.tracked_variables["n_moves"],
        ncmc_mover.N,
    )
)

print(f"Simulating for {equil_md_steps} MD equil steps")
simulation.step(equil_md_steps)

ncmc_mover.reset()  # Reset stats
print("\nProduction....")

print("Running NCMC...")
n_acc = 0
for i in range(2000):
    simulation.step(1000)
    ncmc_mover.move(simulation.context, 1)
    ncmc_mover.report(simulation)
    if ncmc_mover.n_accepted == n_acc + 1:
        n_acc += 1
        moves_since_accepted = 0
    elif ncmc_mover.n_accepted == n_acc:
        moves_since_accepted += 1

    if moves_since_accepted == 350:
        break  # Quit after 350 with no accepted moves


# Setup the output files
# Move ghost waters out of the simulation cell
trj = grand.utils.shift_ghost_molecules(
    ghost_file="ncmc-ghost-ligs.txt",
    topology="MUP1WithGhosts.pdb",
    trajectory="MUP1NCMCholo_raw.dcd",
)

# Recentre the trajectory on a particular residue
# trj = grand.utils.recentre_traj(t=trj, resname='ILE', resid=221)

# Align the trajectory to the protein
grand.utils.align_traj(t=trj, output="MUP1NCMCholo.dcd", reference="MUP1WithGhosts.pdb")

# Write out a trajectory of the GCMC sphere
grand.utils.write_sphere_traj(
    radius=sphere_rad,
    ref_atoms=ref_atoms,
    topology="MUP1WithGhosts.pdb",
    trajectory="MUP1NCMCholo.dcd",
    output="gcmc_sphere.pdb",
    initial_frame=True,
)

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
pdb = PDBFile("../../../../bCD-Equiled.pdb")
lig_pdb = glob.glob("../../*.pdb")
print(lig_pdb)
lig_xml = glob.glob("../../*.xml")
# Add ghost waters,
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(
    pdb.topology, pdb.positions, molfile=lig_pdb[0], n=10, pdb="betaCDWithghosts.pdb"
)
# Load force field and create system
ff = ForceField("../../../../bCD_FINAL.xml", "amber14/tip3p.xml", lig_xml[0])
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=12.0 * angstroms,
    switchDistance=10.0 * angstroms,
    constraints=HBonds,
)

# Langevin Integrator
integrator = BAOABIntegrator(298 * kelvin, 1.0 / picosecond, 0.002 * picoseconds)

# Sphere data
ref_atoms = [
    {"name": "C34", "resname": "MGO", "resid": "1"},
    {"name": "C7", "resname": "MGO", "resid": "1"},
]

sphere_rad = 5.0 * angstrom


# Set up NCMC Sphere sampler
print("Setting up sphere sampler object...")
ncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(
    system=system,
    topology=pdb.topology,
    temperature=298 * kelvin,
    integrator=integrator,
    excessChemicalPotential=None,
    adams=-17.54279661956912,
    standardVolume=None,
    referenceAtoms=ref_atoms,
    sphereRadius=sphere_rad,
    nPertSteps=1499,
    nPropStepsPerPert=50,
    resname="L01",
    ghostFile="ncmc-ghost-ligs.txt",
    log="betaCD.log",
    dcd="betaCD_raw.dcd",
    overwrite=True,
    maxN=1,
    recordTraj=False,
)


# Get platform
platform = Platform.getPlatformByName("CUDA")
platform.setPropertyDefaultValue("Precision", "mixed")
print(f"Running on {platform} platform.")
# Set up simulation object
simulation = Simulation(pdb.topology, system, ncmc_mover.integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(298 * kelvin)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

simulation.reporters.append(
    StateDataReporter(
        stdout,
        500,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        density=True,
        volume=True,
    )
)

print("Initialising...")
# Switch off ghosts
ncmc_mover.initialise(simulation.context, simulation, ghosts)
ncmc_mover.deleteMoleculesInGCMCSphere()

print("Minimising...")
simulation.minimizeEnergy(
    tolerance=0.0001 * kilojoules_per_mole / nanometer, maxIterations=10000
)  # Quick Minimisation

print("Equilibration...")
for i in range(150):
    simulation.step(500)
    ncmc_mover.move(simulation.context, 1)

print(
    "{}/{} equilibration NCMC moves accepted. N = {}".format(
        ncmc_mover.tracked_variables["n_accepted"],
        ncmc_mover.tracked_variables["n_moves"],
        ncmc_mover.N,
    )
)

ncmc_mover.reset()  # Reset stats
print("\nProduction....")

print("Running NCMC...")
for i in range(2000):
    simulation.step(500)
    ncmc_mover.move(simulation.context, 1)
    ncmc_mover.report(simulation)


with open("Insertion_Works.txt", "w") as fi:
    for work in ncmc_move.tracked_variables["insert_works"]:
        fi.write(f"{work}\n")

with open("Deletion_Works.txt", "w") as fi:
    for work in ncmc_move.tracked_variables["delete_works"]:
        fi.write(f"{work}\n")


# Setup the output files
# Move ghost waters out of the simulation cell
trj = grand.utils.shift_ghost_molecules(
    ghost_file="ncmc-ghost-ligs.txt",
    topology="betaCDWithghosts.pdb",
    trajectory="betaCD_raw.dcd",
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

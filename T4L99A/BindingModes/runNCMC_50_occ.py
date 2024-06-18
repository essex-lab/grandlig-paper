import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
import grand
from sys import stdout
import glob




# Load in PDB
pdb = PDBFile("../NPT_equil_no_repulsions.pdb")
lig_pdb = glob.glob("../Me*.pdb")
lig_xml = glob.glob("../*.xml")
print(lig_pdb)
# Add ghost waters,
pdb.topology, pdb.positions, ghosts = grand.utils.add_ghosts(
    pdb.topology,
    pdb.positions,
    molfile=lig_pdb[0],
    n=20,
    pdb="T4Withghosts.pdb",
)
# Load force field and create system
ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", lig_xml[0])
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=12.0 * angstroms,
    switchDistance=10.0 * angstroms,
    constraints=HBonds,
)

# Langevin Integrator
integrator = BAOABIntegrator(
    298 * kelvin, 1.0 / picosecond, 0.002 * picoseconds
)

# Sphere data
ref_atoms = [
    {"name": "CG", "resname": "LEU", "resid": "85"},
    {"name": "CB", "resname": "ALA", "resid": "100"},
]
sphere_rad = 8.0 * angstrom

# Set up NCMC Sphere sampler
print("Setting up sphere sampler object...")
ncmc_mover = grand.samplers.NonequilibriumGCMCSphereSampler(
    system=system,
    topology=pdb.topology,
    temperature=298 * kelvin,
    integrator=integrator,
    excessChemicalPotential=None,
    adams=-7.340809557310219,
    standardVolume=None,
    referenceAtoms=ref_atoms,
    sphereRadius=sphere_rad,
    nPertSteps=499,
    nPropStepsPerPert=50,
    resname="L02",
    ghostFile="ncmc-ghost-ligs.txt",
    log="T4NCMC.log",
    dcd="T4NCMCholo_raw.dcd",
    overwrite=True,
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

# Switch off ghosts
ncmc_mover.initialise(simulation.context, ghosts)
ncmc_mover.deleteMoleculesInGCMCSphere()

print("Minimising...")
simulation.minimizeEnergy(
    tolerance=0.0001 * kilojoule / mole, maxIterations=100000
)  # Quick Minimisation

print("Equilibration...")
for i in range(150):
    simulation.step(500)
    ncmc_mover.move(simulation.context, 1)

print(
    "{}/{} equilibration NCMC moves accepted. N = {}".format(
        ncmc_mover.n_accepted, ncmc_mover.n_moves, ncmc_mover.N
    )
)

ncmc_mover.reset()  # Reset stats
print("\nProduction....")

print("Running NCMC...")
for i in range(10000):
    ncmc_mover.deleteMoleculesInGCMCSphere()
    simulation.step(500)
    ncmc_mover.move(simulation.context, 1)
    simulation.step(500)
    ncmc_mover.report(simulation)

# Setup the output files
# Move ghost waters out of the simulation cell
trj = grand.utils.shift_ghost_waters(
    ghost_file="ncmc-ghost-ligs.txt",
    topology="T4Withghosts.pdb",
    trajectory="T4NCMCholo_raw.dcd",
)

# Recentre the trajectory on a particular residue
# trj = grand.utils.recentre_traj(t=trj, resname='ILE', resid=221)

# Align the trajectory to the protein
grand.utils.align_traj(
    t=trj, output="T4NCMCholo.dcd", reference="T4Withghosts.pdb"
)

# Write out a trajectory of the GCMC sphere
grand.utils.write_sphere_traj(
    radius=sphere_rad,
    ref_atoms=ref_atoms,
    topology="T4Withghosts.pdb",
    trajectory="T4NCMCholo.dcd",
    output="gcmc_sphere.pdb",
    initial_frame=True,
)

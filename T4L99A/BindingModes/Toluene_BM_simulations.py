import numpy as np
import argparse
from openmm.app import *
from openmm import *
from openmm.unit import *
from openmmtools.integrators import BAOABIntegrator
import grandlig as grand
from sys import stdout
import glob


# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("--gpuid", help="Input gpu ID to run on", type=int, default=0)
args = parser.parse_args()

# Set MD and GCNCMC variables
nonbondedMethod = PME
nonbondedCutoff = 1.2 * nanometers
switchDistance = 10.0 * angstroms
ewaldErrorTolerance = 0.0005  # Default
constraints = HBonds
rigidWater = True  # Default

# Use HMR to speed things up
hydrogenMass = 2 * amu

# Integration Options
dt = 0.004 * picoseconds
temperature = 298 * kelvin
friction = 1.0 / picosecond

# GCNCMC params
mu_ex = None
concentration = None
volume_per_lig = None

switching_time = 50 * picoseconds

# Because we're using HMR we can halve n_prop for the same amount of time between lambdas
# Alternatively, add more MD between and use less lambdas (quicker, but not as good.)
# NOTE: Future studies will look at optimising the switching times

n_prop = 25
n_pert = int((switching_time / (n_prop * dt)) - 1)

adams_value = (
    -7.340809557310219
)  # Explicitly set adams value to B50 as reported in paper

# Load in PDB
pdb = PDBFile("../T4NoCoSolv.pdb")
lig_pdb = glob.glob("../Me*.pdb")
lig_xml = glob.glob("../*.xml")

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
    nonbondedMethod=nonbondedMethod,
    nonbondedCutoff=nonbondedCutoff,
    switchDistance=switchDistance,
    constraints=constraints,
    hydrogenMass=hydrogenMass,
    rigidWater=rigidWater,
)

# Langevin Integrator
integrator = BAOABIntegrator(temperature, friction, dt)

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
    excessChemicalPotential=mu_ex,
    adams=adams_value,
    standardVolume=volume_per_lig,
    referenceAtoms=ref_atoms,
    sphereRadius=sphere_rad,
    nPertSteps=n_pert,
    nPropStepsPerPert=n_prop,
    resname="L02",
    ghostFile="ncmc-ghost-ligs.txt",
    log="T4NCMC.log",
    dcd="T4NCMCholo_raw.dcd",
    overwrite=True,
    maxN=1,  # NOTE: Setting maxN=1 will instantly reject any attempted insertion move in order to save time. Only use when certain your site fits one ligand
)

# Get platform
platform = Platform.getPlatformByName("CUDA")
platform.setPropertyDefaultValue("Precision", "mixed")
platform.setPropertyDefaultValue("DeviceIndex", f"{args.gpuid}")
print(
    f"Running on the {platform.getName()} platform. On device with index: {platform.getPropertyDefaultValue('DeviceIndex')}"
)

# Set up simulation object
simulation = Simulation(pdb.topology, system, ncmc_mover.integrator, platform)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature)
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
ncmc_mover.initialise(simulation.context, simulation, ghosts)
ncmc_mover.deleteMoleculesInGCMCSphere()

print("Minimising...")
simulation.minimizeEnergy(maxIterations=100000)  # Quick Minimisation

print("Equilibration...")
for i in range(100):
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
for i in range(10000):
    simulation.step(500)
    ncmc_mover.move(simulation.context, 1)
    ncmc_mover.report(simulation, data=True)

# Setup the output files
# Move ghost waters out of the simulation cell
trj = grand.utils.shift_ghost_molecules(
    ghost_file="ncmc-ghost-ligs.txt",
    topology="T4Withghosts.pdb",
    trajectory="T4NCMCholo_raw.dcd",
)


# Align the trajectory to the protein
grand.utils.align_traj(t=trj, output="T4NCMCholo.dcd", reference="T4Withghosts.pdb")

# Write out a trajectory of the GCMC sphere
grand.utils.write_sphere_traj(
    radius=sphere_rad,
    ref_atoms=ref_atoms,
    topology="T4Withghosts.pdb",
    trajectory="T4NCMCholo.dcd",
    output="gcmc_sphere.pdb",
    initial_frame=True,
)

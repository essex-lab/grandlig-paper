import math

import numpy as np
import argparse
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
import grand
from sys import stdout
import glob


def read_moves(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    n_completed_lig = int(lines[-1].split()[4])
    n_accepted_lig = int(lines[-1].split()[7].strip("("))

    n_completed_wat = int(lines[-2].split()[4])
    n_accepted_wat = int(lines[-2].split()[7].strip("("))

    return n_completed_lig, n_accepted_lig, n_completed_wat, n_accepted_wat


# Check which run this is
parser = argparse.ArgumentParser()
parser.add_argument(
    "-r",
    "--run",
    type=int,
    default=2,
    help="Which leg this represents in the full simulation",
)
args = parser.parse_args()

# Input Files
pdb = PDBFile(
    "sol-ghosts2.pdb"
)  # Because this is a continuation, load in the file with ghosts - only needed for topology
lig_pdb = glob.glob("../*Ace*.pdb")
lig_xml = glob.glob("../*Ace*.xml")

# System Configuration
nonbondedMethod = PME
nonbondedCutoff = 1.2 * nanometers
switchDistance = 10.0 * angstroms
ewaldErrorTolerance = 0.0005  # Default
constraints = HBonds
rigidWater = True  # Default
constraintTolerance = 0.000001
hydrogenMass = 1.0 * amu

# Integration Options
dt = 0.002 * picoseconds
temperature = 298 * kelvin
friction = 1.0 / picosecond
"""
Dont need these for GCMC
pressure = 1.0*atmospheres
barostatInterval = 25
"""

# GCNCMC params and variables
mu_ex = -3.25 * kilocalories_per_mole
volume_per_lig = 3360 * angstroms**3
n_pert = 499
n_prop = 50

wat_mu_ex = -6.09 * kilocalories_per_mole
wat_vol = 31.4788 * angstroms**3
n_pert_wat = 99
n_prop_wat = 50

# Loading some old variables
if args.run == 2:
    wat_ghosts = grand.utils.read_ghosts_from_file("hoh-ghosts.txt")[-1]
    n_moves_lig, n_accepted_lig, n_moves_hoh, n_accepted_hoh = read_moves("hoh.log")
    print(n_moves_lig, n_accepted_lig, n_moves_hoh, n_accepted_hoh)
    ghosts = grand.utils.read_ghosts_from_file("ncmc-ghost-ligs.txt")[-1]
else:
    n_moves_lig, n_accepted_lig, n_moves_hoh, n_accepted_hoh = read_moves(
        "hoh-{}.log".format(args.run - 1)
    )
    wat_ghosts = grand.utils.read_ghosts_from_file(
        "hoh-ghosts-{}.txt".format(args.run - 1)
    )[-1]
    ghosts = grand.utils.read_ghosts_from_file(
        "ncmc-ghost-ligs-{}.txt".format(args.run - 1)
    )[-1]


# Load in the .rst7 to get the checkpointed positions and velocities
if args.run == 2:
    rst7 = AmberInpcrdFile("solution.rst7")
else:
    rst7 = AmberInpcrdFile("solution-{}.rst7".format(args.run - 1))


# Get platform
platform = Platform.getPlatformByName("CUDA")
platformProperties = {"Precision": "mixed"}
platform.setPropertyDefaultValue("Precision", "mixed")
print(f"Running on {platform} platform.")


# Load force field and create system
ff = ForceField("amber14-all.xml", "amber14/tip3p.xml", lig_xml[0])
system = ff.createSystem(
    pdb.topology,
    nonbondedMethod=nonbondedMethod,
    nonbondedCutoff=nonbondedCutoff,
    switchDistance=switchDistance,
    constraints=constraints,
    constraintTolerance=constraintTolerance,
    rigidWater=rigidWater,
    hydrogenMass=hydrogenMass,
)

# Langevin Integrator
integrator = BAOABIntegrator(temperature, friction, 0.002 * picoseconds)

# Set up NCMC System samplers
# Set up the forces together?
params, custom_nb_force, ele_except, vdw_except = grand.utils.create_custom_forces(
    system, pdb.topology, ["L02", "HOH"]
)


###### WATER #######
# Create GCMC sampler for water
water_sampler = grand.samplers.NonequilibriumGCMCSystemSampler(
    system=system,
    topology=pdb.topology,
    temperature=298 * kelvin,
    integrator=integrator,
    nPertSteps=n_pert_wat,
    nPropStepsPerPert=n_prop_wat,
    timeStep=dt,
    # Make sure the two parameters below are correct for the ligand concentration
    excessChemicalPotential=wat_mu_ex,
    standardVolume=wat_vol,
    boxVectors=np.array(pdb.topology.getPeriodicBoxVectors()),
    createCustomForces=False,
    log="hoh-{}.log".format(args.run),
    rst="solution-{}.rst7".format(args.run),
    ghostFile="hoh-ghosts-{}.txt".format(args.run),
    overwrite=True,
)

water_sampler.setCustomForces(params["HOH"], custom_nb_force, ele_except, vdw_except)

print("Setting up system Ligand sampler object...")
ncmc_mover = grand.samplers.NonequilibriumGCMCSystemSampler(
    system=system,
    topology=pdb.topology,
    temperature=temperature,
    integrator=integrator,
    excessChemicalPotential=mu_ex,
    adams=None,
    standardVolume=volume_per_lig,
    nPertSteps=n_pert,
    nPropStepsPerPert=n_prop,
    timeStep=dt,
    boxVectors=np.array(pdb.topology.getPeriodicBoxVectors()),
    resname="L02",
    createCustomForces=False,
    ghostFile="ncmc-ghost-ligs-{}.txt".format(args.run),
    log="bulk_conc_0.1-{}.log".format(args.run),
    rst="checkpoint-{}.rst7".format(args.run),
    dcd="bulkconc_raw-{}.dcd".format(args.run),
    overwrite=True,
)
ncmc_mover.setCustomForces(params["L02"], custom_nb_force, ele_except, vdw_except)

# Set up simulation object
simulation = Simulation(
    pdb.topology,
    system,
    ncmc_mover.integrator,
    platform,
    platformProperties=platformProperties,
)
simulation.context.setPositions(rst7.getPositions())  # From Checkpoint
simulation.context.setVelocitiesToTemperature(rst7.getVelocities())  # From checkpoint
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
water_sampler.initialise(simulation.context, simulation, wat_ghosts)


# Set the number of moves to that left off at
water_sampler.n_moves = n_moves_hoh
water_sampler.n_accepted = n_accepted_hoh
ncmc_mover.n_moves = n_moves_lig
ncmc_mover.n_accepted = n_accepted_lig


print("Minimising...")
simulation.minimizeEnergy(
    tolerance=0.0001 * kilojoule / mole, maxIterations=10000
)  # Quick Minimisation

print("Equilibration...")
# Run some MD for equil
equil_md_steps = math.ceil((0.2 * nanosecond) / (0.002 * picoseconds))
print(f"Simulating for {equil_md_steps} MD equil steps at 2 fs")
simulation.step(equil_md_steps)


print("\nProduction....")

print("Running NCMC...")
while ncmc_mover.n_moves < 100000:
    ncmc_mover.move(simulation.context, 1)
    water_sampler.move(simulation.context, 3)
    simulation.step(10000)
    water_sampler.report(simulation)
    ncmc_mover.report(simulation)

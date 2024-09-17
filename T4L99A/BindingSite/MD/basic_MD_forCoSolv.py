from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from openmmtools.integrators import BAOABIntegrator
import argparse
from mdtraj.reporters import DCDReporter         # <-- new import from mdtraj


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--pdb', type=str)
parser.add_argument('-x', '--xml', type=str)
parser.add_argument('--hmr', help='Use HMR?', action='store_true')
args = parser.parse_args()



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


pdb = PDBFile(args.pdb)

# Create system
forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml', args.xml)
system = forcefield.createSystem(pdb.topology,
                                 nonbondedMethod=nonbondedMethod,
                                 nonbondedCutoff=nonbondedCutoff,
                                 switchDistance=switchDistance,
                                 rigidWater=rigidWater,
                                 constraints=constraints,
                                 hydrogenMass=hydrogenMass)


barostat = MonteCarloBarostat(1*bar, temperature, 25)

# Define integrator
integrator = BAOABIntegrator(temperature, friction, dt)

# Define platform and set precision
platform = Platform.getPlatformByName('CUDA')
platform.setPropertyDefaultValue('Precision', 'mixed')

# Create simulation object
simulation = Simulation(pdb.topology, system, integrator, platform)

# Set positions, velocities and box vectors
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperature)
simulation.context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())
simulation.reporters.append(StateDataReporter('simulation_log.txt', 1000, step=True, time=True,
        potentialEnergy=True, temperature=True, density=True, volume=True, elapsedTime=True))


# Minimise the system 
print('Minimising...')
simulation.minimizeEnergy()
print('Minimised.')

# Equilibrate the system for 1 ns NVT
print('Equilibrating NVT...')
n_steps = (1 * nanosecond) / dt
simulation.step(int(n_steps))
print('NVT Equilibrated.')


# Add barostat and minimise more
system.addForce(barostat)
simulation.context.reinitialize(preserveState=True)
print('Equilibrating NPT...')
n_steps = (4 * nanosecond) / dt
simulation.step(int(n_steps))

print('NPT Equilibrated.')


print("Production...")
simulation.reporters.append(DCDReporter('cosolv_traj.dcd', 2500))

n_steps = (40 * nanosecond) / dt
simulation.step(int(n_steps))

simulation.reporters.append(DCDReporter('final_10ns.dcd', 2500))
n_steps = (10 * nanosecond) / dt
simulation.step(int(n_steps))

print("Production Complete")


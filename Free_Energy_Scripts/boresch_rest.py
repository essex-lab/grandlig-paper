from openmm.unit import *
import numpy as np
import mdtraj as md
from openmmtools.constants import STANDARD_STATE_VOLUME
from openmm import CustomCompoundBondForce

def BoreschRestraint(pdb, restrained_atoms=None, K_r=None, K_theta=None, K_phi=None, r=None, a=None, d=None,
                     lambda_controlled='lambda_boresch', kT=None):
    """
    Function to create a CustomCompoundBond OpenMM force for Boresch restraints.
    Parameters
    ----------
    pdb : openn.app.PDBFile
    restrained_atoms : list of ints
        list of atoms to restrain in for [l3 l2 l1 p1 p2 p3]
    K_r :
    K_theta
    K_phi
    r : float
        reference distance to use
    a : list of floats
        Reference angles to use
    d : list of floats
        Reference dihedral values to use
    lambda_controlled
    kT : openmm.unit.Quantity
        Value of kT for the simulation

    Returns
    -------

    """
    print('Adding Boresch restraints..')
    if K_r is None:
        K_r = 10 * kilocalories_per_mole / angstroms**2
    if K_theta is None:
        K_theta = 10 * kilocalories_per_mole / radians**2
    if K_phi is None:
        K_phi = 10 * kilocalories_per_mole / radians**2

    energy_function = 'E;' \
                      'E = (step(distance(p3,p4) - r_aA0) * (K_r/2)*(distance(p3,p4) - r_aA0)^2)' \
                      '+ (K_thetaA/2)*(angle(p2,p3,p4)-theta_A0)^2 + (K_thetaB/2)*(angle(p3,p4,p5)-theta_B0)^2' \
                      '+ (K_phiA/2)*dphi_A^2 + (K_phiB/2)*dphi_B^2 + (K_phiC/2)*dphi_C^2; ' \
                      'dphi_A = dA - floor(dA/(2*pi)+0.5)*(2*pi); dA = dihedral(p1,p2,p3,p4) - phi_A0; ' \
                      'dphi_B = dB - floor(dB/(2*pi)+0.5)*(2*pi); dB = dihedral(p2,p3,p4,p5) - phi_B0; ' \
                      'dphi_C = dC - floor(dC/(2*pi)+0.5)*(2*pi); dC = dihedral(p3,p4,p5,p6) - phi_C0; ' \
                      'pi = {0}; '.format(np.pi)

    if lambda_controlled:
        energy_function = f'{lambda_controlled} * ' + energy_function

    positions = pdb.getPositions(asNumpy=True)
    topology = pdb.topology
    pdb_atoms = [atom for atom in pdb.topology.atoms()]
    t = md.Trajectory(positions, topology)  # Define mdtraj trajectory
    # Distance
    atom_pairs = [restrained_atoms[2:4]]
    if r is None:
        r_aA0 = md.geometry.compute_distances(t, atom_pairs, periodic=False)[0][0] * nanometers
    else:
        r_aA0 = r * nanometers
    print(f'raA0 = {r_aA0} between atoms {pdb_atoms[atom_pairs[0][0]]} and {pdb_atoms[atom_pairs[0][1]]}')

    # Angles
    atom_triplets = [restrained_atoms[i:(i + 3)] for i in range(1, 3)]
    if a is None:
        angles = md.geometry.compute_angles(t, atom_triplets, periodic=False)[0]
        theta_A0 = angles[0] * radians
        theta_B0 = angles[1] * radians
    else:
        theta_A0 = (a[0] * degrees).in_units_of(radians)
        theta_B0 = (a[1] * degrees).in_units_of(radians)
    print(f'theta_A0 = {theta_A0.in_units_of(degrees)} between {[pdb_atoms[i] for i in atom_triplets[0]]}')
    print(f'theta_B0 = {theta_B0.in_units_of(degrees)} between {[pdb_atoms[i] for i in atom_triplets[1]]}')

    # Dihedrals
    atom_quadruplets = [restrained_atoms[i:(i + 4)] for i in range(3)]
    if d is None:
        dihedrals = md.geometry.compute_dihedrals(t, atom_quadruplets, periodic=False)[0]
        phi_A0 = dihedrals[0] * radians
        phi_B0 = dihedrals[1] * radians
        phi_C0 = dihedrals[2] * radians
    else:
        phi_A0 = (d[0] * degrees).in_units_of(radians)
        phi_B0 = (d[1] * degrees).in_units_of(radians)
        phi_C0 = (d[2] * degrees).in_units_of(radians)

    print(f'phi_A0 = {phi_A0} between {[pdb_atoms[i] for i in atom_quadruplets[0]]}')
    print(f'phi_B0 = {phi_B0} between {[pdb_atoms[i] for i in atom_quadruplets[1]]}')
    print(f'phi_C0 = {phi_C0} between {[pdb_atoms[i] for i in atom_quadruplets[2]]}')

    parameter_items = {'K_r': K_r, 'K_thetaA': K_theta, 'K_thetaB': K_theta,
                       'K_phiA': K_phi, 'K_phiB': K_phi, 'K_phiC': K_phi,
                       'r_aA0': r_aA0, 'theta_A0': theta_A0,
                       'theta_B0': theta_B0, 'phi_A0': phi_A0, 'phi_B0': phi_B0,
                       'phi_C0': phi_C0}  # Dictionary of the parameters and values

    for parameter in parameter_items.keys():
        print(f'{parameter} = {parameter_items[parameter]}')  # Print to the user
        # Add the values to the energy function as these are parameters which dont change throughout the simulation
        energy_function += f'{parameter}={parameter_items[parameter].value_in_unit_system(md_unit_system)}; '

    # Setup and create the force
    n_particles = 6
    restraint_force = CustomCompoundBondForce(n_particles, energy_function)
    if lambda_controlled:
        restraint_force.addGlobalParameter(lambda_controlled, 1.0)
    restraint_force.addBond(restrained_atoms, [])
    restraint_force.setUsesPeriodicBoundaryConditions(True)

    # Calculate the analytical correction for realising the restraints
    t_A0rads = theta_A0.in_units_of(radians)._value
    t_B0rads = theta_B0.in_units_of(radians)._value
    numerator = 8.0 * np.pi**2 * STANDARD_STATE_VOLUME * np.sqrt(K_r * (K_theta**2) * (K_phi**3)) * 1*radians**5
    denominator = r_aA0**2 * np.sin(t_A0rads) * np.sin(t_B0rads) * (2 * np.pi * kT)**3
    std_state_correction = -kT * np.log(numerator/denominator)
    return restraint_force, std_state_correction

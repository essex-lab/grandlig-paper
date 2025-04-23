#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to represent the positions of the ligands onto a grid to visually predict free energy
Takes trajectory of the simulation and calculates the mean position of each ligand in each frame
Then finds out which tick mark on the axis the mean is closest two and then adds +1 to that position to
show that a ligand was there.

Writes out the data into a .cube file that can be read by visualisation software to generate an isosurface
- Will Poole
"""


import argparse
import numpy as np
from tqdm import tqdm
import mdtraj
from numba import jit, cuda
from simtk import unit

#@staticmethod@
@jit(target_backend="cuda")
def get_frame_grid(lig_atoms, threshold, grid_shape, xyz):
    """
    Function which takes a cube file and outputs the essential information lines to a new cube file (first ~7 lines)

    Parameters
    ----------
    lig_atoms : str
        List containing all the atoms from all the ligands with a corresponding atom number
    threshold : float
        The given distance threshold for a probe to be within to a grid point to be counted to that grid point
    grid_shape : np.shape
        The shape of the main grid containing the data
    xyz : float (np array)
        Coordinates of the ligand atoms at a certain frame (given by the loops
    Returns
    -------
    n : np array
       Returns the grid that has been calculated for this frame in the loop
    """
    # Get the grid occupancy for a single frame
    # Count this frame only
    grid_frame = np.zeros(grid_shape)  # Creates a grid specific to the frame in the loop based on main grid
    # Loop over ligand atoms
    for index in lig_atoms:  # For each atom number in the list of ligand atoms
        for i, x in enumerate(xs):
            dx = abs(x - 10*xyz[index, 0])
            if dx > threshold:
                continue
            for j, y in enumerate(ys):
                dy = abs(y - 10*xyz[index, 1])
                if dy > threshold:
                    continue
                for k, z in enumerate(zs):
                    dz = abs(z - 10*xyz[index, 2])
                    if dz < threshold:
                        if np.linalg.norm(np.array([dx, dy, dz])) < threshold:
                            grid_frame[i, j, k] = 1
    return grid_frame


# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-di", "--dcdin", help="Input .dcd Simulation File", default='system.dcd')
parser.add_argument("-p", "--pdb", help="Input PDB file", default='system.pdb')
parser.add_argument("-r", "--resname", help="Provide a the resname of the ligand")
parser.add_argument("-y", "--origin", help="Do you want to supply your own origin details?",
                    default=None)
parser.add_argument("-t", "--time", help="Simulation time to start the grid analysis.", type=int, default=0)
parser.add_argument("-b", "--bulk", help="Bulk occupancy.", type=float, default=None)
parser.add_argument("-o", "--output", help="Provide a filename for the output grid.", default='grid.cube')
parser.add_argument("-dt", "--dt", help="Enter the time between each frame in ps.", default=20, type=float)


args = parser.parse_args()

# Creating the MDTraj "Universe"
t = mdtraj.load(args.dcdin, top=args.pdb)  # Loads trajectory
n_frames, n_atoms, n_dims = t.xyz.shape  # Extracts useful information
tot_time = np.round((n_frames * args.dt / 1000)) # frames * time between each frame (ps) / 1000 to get into ns
print(f'Total simulation time: {tot_time} ns in {n_frames} frames.')
start_time_wanted = args.time
frames_per_ns = n_frames / tot_time
start_frame = start_time_wanted * frames_per_ns
n_frames_new = n_frames - start_frame
print(f'Using {n_frames_new} frames for the analysis for total of {(n_frames_new * args.dt) / 1000} ns of data starting from {args.time} ns.')

b = 0.529177
if args.origin:  # If you supply origin details this bit of code will extract those details - useful for building grids
    min_coords = [0, 0, 0]
    max_coords = [0, 0, 0]
    with open(args.origin, 'r') as f1:
        lines = f1.readlines()[2:7]
        for i in range(3):
            min_coords[i] = float(lines[0].split()[i + 1]) * b
            max_coords[i] = (float(lines[0].split()[i + 1]) +
                             ((float(lines[i + 1].split()[0])) * float(lines[i+1].split()[i+1]))) * b
        x_tick, y_tick, z_tick = int(lines[1].split()[0]), int(lines[2].split()[0]), int(lines[3].split()[0])
        x_tickdist, y_tickdist, z_tickdist = (float(lines[1].split()[1]) * b), (float(lines[2].split()[2]) * b), \
                                             (float(lines[3].split()[3]) * b)
    xs = np.linspace(min_coords[0], max_coords[0], int(x_tick))
    ys = np.linspace(min_coords[1], max_coords[1], int(y_tick))
    zs = np.linspace(min_coords[2], max_coords[2], int(z_tick))

else:  # If NOT supplied then we must extract that data ourselves.
    protein_ids = t.topology.select('protein')  # Selects protein atoms and store the atom number in an array
    # Working out the maximum and minimum coordinates that the protein has in order to define the grid size
    max_coords = [-1e8, -1e8, -1e8]   # Creates a list with 3 entries of a v.small number
    min_coords = [1e8, 1e8, 1e8]  # Same but for a v.large number

    for index in tqdm(protein_ids):  # Loops over the array of protein atom numbers
        for i in range(3):  # Loops over the x, y, z dimensions
            max_tmp = np.max(t.xyz[:, index, i])   # Extracts the max coordinate in x y and z for specific atom
            min_tmp = np.min(t.xyz[:, index, i])  # Same for min, across all frames
            if max_tmp > max_coords[i]:  # If the temp max is greater than what is already in the list then
                max_coords[i] = max_tmp  # Change the max coord to the max temp
            if min_tmp < min_coords[i]:  # if temp min is less than what is already in list
                min_coords[i] = min_tmp  # Change min to min temp

    for i in range(3):  # Loops over each dimension and multiplies the min/max value by 10 to go from nm to A and then
                        # plus and minus 2 from max and min to give a 2A buffer
        min_coords[i] *= 10
        max_coords[i] *= 10
        min_coords[i] -= 2
        max_coords[i] += 2

    # Works out the difference between the x, y, z axis limits of the grid - Divide this by how many tick marks
    # you would want along the axis - In this case we want ~3A per tick so /1.
    x_tick = 3 * (int((max_coords[0])-(min_coords[0])))
    y_tick = 3 * (int((max_coords[1])-(min_coords[1])))
    z_tick = 3 * (int((max_coords[2])-(min_coords[2])))

    # Creates the axis values
    xs = np.linspace(min_coords[0], max_coords[0], x_tick)
    ys = np.linspace(min_coords[1], max_coords[1], y_tick)
    zs = np.linspace(min_coords[2], max_coords[2], z_tick)
    

    # Calculates the distance between the tick marks (axis values)
    x_tickdist = xs[1] - xs[0]
    y_tickdist = ys[1] - ys[0]
    z_tickdist = zs[1] - zs[0]

print(min_coords, max_coords)


grid = np.zeros((x_tick, y_tick, z_tick))  # Creates numpy matrix of zeros using the grid tick marks

threshold = 1.6  # If something is within this distance, it matches a point

# Selects the atoms in the ligands
lig_atoms = t.topology.select("resn {} and not name 'H*'".format(args.resname))

# Loop over frames
for f in tqdm(range(n_frames)):
    # Add onto main grid
    grid += get_frame_grid(lig_atoms, threshold, grid.shape, t.xyz[f, :, :])
grid /= n_frames  # Divides every value in the grid by the number of frames to normalise the data
RT = 1.9872042e-3 * 300 # in kcals


if args.bulk:
    grid = -(RT * np.log(grid/args.bulk))
    grid[grid == np.inf] = 0 # change the infs back to 0
    grid[grid == -np.inf] = 0


with open(args.output, 'w') as out:  # Opens up the grid file as full write
    out.write(" CPMD CUBE FILE\n OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z\n")  # First and second line = Comments
    out.write("{:>5d}{:>12f}{:>12f}{:>12f}\n".format(1, (min_coords[0]/b),
                                                     (min_coords[1])/b,
                                                     (min_coords[2])/b))  # Third line = Origin
    out.write("{:>5d}{:>12f}{:>12f}{:>12f}\n".format(x_tick, (x_tickdist/b), 0.000000, 0.000000))
    # Fourth line = xtick dists
    out.write("{:>5d}{:>12f}{:>12f}{:>12f}\n".format(y_tick, 0.000000, y_tickdist/b, 0.000000))
    out.write("{:>5d}{:>12f}{:>12f}{:>12f}\n".format(z_tick, 0.000000, 0.000000, z_tickdist/b))
    out.write("{:>5d}{:>12f}{:>12f}{:>12f}{:>12f}\n".format(1, 1.000000, 0.000000, 0.000000, 0.000000))
    # For loops to write out the data in the formatting required for a cube file - See online tutorial
    for x in tqdm(range(len(xs))):
        for y in range(len(ys)):
            for z in range(len(zs)):
                out.write("  {:>2.5E}".format(grid[x, y, z],))
                if z % 6 == 5:
                    out.write("\n")

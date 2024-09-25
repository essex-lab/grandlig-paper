
import argparse
import numpy as np
from matplotlib import pyplot as plt
from simtk.unit import *


def calc_concentration(N, V):
    """
    Calculate the concentration for a given volume and number of molecules

    Parameters
    ----------
    N : int
        Number of molecules
    V : simtk.unit.Quantity
        Volume, in appropriate units

    Returns
    -------
    conc : simtk.unit.Quantity
        Calculated concentration, in units of M
    """
    n_moles = N / AVOGADRO_CONSTANT_NA
    conc = n_moles / V
    return conc.in_units_of(molar)


parser = argparse.ArgumentParser()
parser.add_argument('-d', '--data', nargs='+', default=['gcmd.log'], help='GCMD log file.')
parser.add_argument('-n', '--nmols', default=None, type=int, help='Initial number of waters in the simulation')
parser.add_argument('-v', '--volume', default=None, type=float, help='Volume of the system (nm^3)')
parser.add_argument('-o', '--output', default='uvt.csv', help='Name of the CSV file containing the output')
parser.add_argument('-r', '--n_repeats', default=1, type=int, help='Number of repeats going into the thing')

args = parser.parse_args()

# Print the sensitivity of the concentration to the number of molecules
print('Sensitivity = {} M'.format(calc_concentration(1, args.volume * nanometers ** 3)._value))

moves = []
move_dict = {}  # Store the concentrations for a given number of moves (proxy for time)

# Calculate the initial concentration
if args.nmols is not None:
   moves.append(0)
   conc = calc_concentration(args.nmols, args.volume * nanometers ** 3)._value
   move_dict[0] = [conc]

# Read in the concentration at each point for each simulation
for filename in args.data:
    with open(filename, 'r') as f:
        for line in f.readlines()[1::2]:
            print(line)
            # Check that this is a reporting line
            if 'Current N' in line:
                for i in range(len(line.split())):
                    # Read the number of moves executed at this point
                    if line.split()[i] == 'move(s)':
                        n_moves = int(line.split()[i-1])
                        if n_moves not in moves:
                            moves.append(n_moves)
                    # Read the value of N at this point
                    elif line.split()[i] == 'Current':
                        n_wats = int(line.split()[i+3].strip('.'))
                # Calculate concentration
                conc = calc_concentration(n_wats, args.volume * nanometers ** 3)._value
                # Add this concentration to the appropriate list (or create a new one)
                try:
                    move_dict[n_moves].append(conc)
                except:
                    move_dict[n_moves] = [conc]

# Ensure the dicts are same length
len_dict = 0
for key in move_dict.keys():
    if len(move_dict[key]) > len_dict:
        len_dict = len(move_dict[key])

for key in move_dict.keys():
    to_add = len_dict - len(move_dict[key])
    for i in range(to_add):
        move_dict[key].append(np.NaN)

import pandas as pd
# col_names = ["move"]
# n_repeats = 8
# for i in range(n_repeats):
#     col_names.append(f"repeat_{i+1}")
df = pd.DataFrame.from_dict(move_dict)
df = df.transpose()

df.index.name = "Moves"
header = [f"repeat_{i}" for i in range(len(df.columns))]
df.to_csv(f"{args.output}.csv", index=True, header=header)
quit()

# Write data to a CSV file
with open(args.output, 'w') as f:
    f.write('Moves,Concentration (M)\n')
    for i in range(len(moves)):
        # Convert this number of moves to a point in time
        f.write('{}'.format(moves[i]))
        for c in move_dict[moves[i]]:
            f.write(',{}'.format(c))
        f.write('\n')


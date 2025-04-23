# Modules to import
import argparse
from simtk.unit import *
import numpy as np
from glob import glob


# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-f", '--folders', help='Input folder names of the fragment repeats', type=str, nargs='+')
parser.add_argument("-fn", '--name', help='Input file name of the log file', type=str)
parser.add_argument("-frag", '--frag', help='Input name of the fragment', type=str, default='Ligand')

args = parser.parse_args()



def get_avgN(file):
    with open(file, 'r') as fi:
        lines = fi.readlines()
        final_line = lines[-1].split()
        avg_N = float(final_line[-1])
        return avg_N


repeats = args.folders
n_repeats = len(repeats)
log_name = args.name

# Get all the B values simulated for a particular ligand. N.B these dont have to be simulated in eavery repeat
unique_Bs = []
for repeat in repeats:
    folders = glob(f'{repeat}/*.*/')
    print(folders)
    for f in folders:
        b = float(f.split('/')[-2])
        if b not in unique_Bs:
            unique_Bs.append(b)

Bs = unique_Bs.copy()
print(Bs)

B_ids = []

for i in range(len(Bs)):
    B_ids.append(i)
Bs.sort()


avg_N_mols = np.zeros((n_repeats, len(B_ids)))  # Dict of B values and the current N waters at each frame
for repeat_id, repeat in enumerate(repeats):
    for i, B in enumerate(Bs):
        print(f'{repeat}/{B}/{log_name}')
        try:
            avg_N = get_avgN(f'{repeat}/{B}/{log_name}')
            avg_N_mols[repeat_id][i] = avg_N
        except FileNotFoundError:
            avg_N_mols[repeat_id][i] = np.NAN


# Write out the data to a csv file!
with open(f'{args.frag}_NOpost_procress_data.csv', 'w') as fi:
    for i in range(len(Bs)):
        fi.write(f'{Bs[i]}')
        for repeat_id in range(n_repeats):
            fi.write(f', {avg_N_mols[repeat_id, i]}')
        fi.write('\n')



import argparse
import MDAnalysis as mda
import numpy as np
import re
from simtk.openmm.app import pdbxfile
import os
from tqdm import tqdm
# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-t", "--traj", help="Input list of trajectory simulation files", default=[], nargs='+')
parser.add_argument("-p", "--pdb", help="Input PDB file", default='system.pdb')
parser.add_argument("-s", "--skip", help="n_frames to skip", default=1, type=int)

args = parser.parse_args()

for t in args.traj:
    print(t)
    print(f'File {t} : {os.path.exists(t)}')



ref = mda.Universe(args.pdb)


all = ref.select_atoms('all')

with mda.Writer('combined-new.dcd', all.n_atoms) as w:
    for traj in args.traj:
        u = mda.Universe(args.pdb, traj, dt=1)
        all_ = u.select_atoms('all')
        for ts in tqdm(u.trajectory[1500::]):
            w.write(all_)


import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from matplotlib import pyplot as plt
import numpy as np
import pickle as pkl
import argparse
from tqdm import tqdm

parser = argparse.ArgumentParser()
parser.add_argument("--pdb")
parser.add_argument("--traj")
parser.add_argument("--resid")
parser.add_argument("--out")

args = parser.parse_args()

u = mda.Universe(args.pdb, args.traj)
n_frames = len(u.trajectory)


val = u.select_atoms(f"resname VAL and resid {args.resid}")
val_res = val.residues[0]
chi1 = [val_res.chi1_selection()]
dihs = dihedrals.Dihedral(chi1).run(verbose=True)


with open(f"{args.out}.txt", "w") as fo:
    for i in tqdm(range(dihs.angles.T.shape[1])):
        fo.write(f"{dihs.angles.T[0][i]}\n")

plt.plot(range(n_frames), dihs.angles.T[0], ls='none', marker='x')
plt.xlabel('Frame')
plt.ylabel('Angle (degrees)')
plt.title("p-xylene / T4L99A - Val111 Dihedral Angle")
plt.savefig(f"{args.out}.pdf", bbox_inches='tight')

plt.clf()
plt.hist(dihs.angles.T[0], density=True)
plt.ylabel('Density')
plt.xlabel('Angle (degrees)')
plt.title("p-xylene / T4L99A - Val111 Dihedral Angle")
plt.savefig(f"{args.out}_hist.pdf", bbox_inches='tight')

#plt.show()
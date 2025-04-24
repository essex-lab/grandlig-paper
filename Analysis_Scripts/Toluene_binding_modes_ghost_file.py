import mdtraj
import MDAnalysis as mda
import numpy as np
import MDAnalysis.core.topologyobjects as mda_topo
from matplotlib import pyplot as plt
import argparse
from tqdm import tqdm

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-p", '--pdb', help='Input file name of the PDB file', type=str)
parser.add_argument("-d", '--dcd', help='Input file name of the DCD file', type=str)

args = parser.parse_args()

small_font = 16
medium_font = 18
large_font = 20

plt.rc('figure', titlesize=large_font)
plt.rc('font', size=small_font)
plt.rc('axes', titlesize=small_font)
plt.rc('axes', labelsize=medium_font)
plt.rc('xtick', labelsize=small_font)
plt.rc('ytick', labelsize=small_font)
plt.rc('legend', fontsize=medium_font)

# Sphere data
ref_atoms = [{'name': 'CG', 'resname': 'LEU', 'resid': '85'},
             {'name': 'CB', 'resname': 'ALA', 'resid': '100'}]
sphere_radius = 8.0


frames_with_lig = []
lig_resid = []

u = mda.Universe(args.pdb, args.dcd)
arg120 = u.select_atoms('resid 120 and name CA')
print(arg120)
#quit()
arg120_id = arg120.indices[0]

#ligs = u.select_atoms('resname L02')

dihedrals = []
for i, frame in tqdm(enumerate(range(len(u.trajectory)))):
    #print(i)
    u.trajectory[frame] # Move there in traj
    lig = u.select_atoms(f'resname L02')
    #print(lig)
    dihedral_lig_atoms = {}
    for atom in lig:
        if atom.name in ['C1', 'C7', 'C5']:
            dihedral_lig_atoms[atom.name] = atom.index
    atoms = [arg120_id, dihedral_lig_atoms['C1'], dihedral_lig_atoms['C7'], dihedral_lig_atoms['C5']]
    #atoms_tester = [dihedral_lig_atoms['C5'], dihedral_lig_atoms['C7'], dihedral_lig_atoms['C1'], arg120_id]
    #print(mda_topo.Dihedral(atoms, u).dihedral(), mda_topo.Dihedral(atoms_tester, u).dihedral())
    dihedrals.append(mda_topo.Dihedral(atoms, u).dihedral())

dihedrals = np.asarray(dihedrals)
#dihedrals += 180
print(dihedrals[0:10])
dihedrals_rad = (dihedrals * np.pi) / 180
print(dihedrals_rad[0:10])
B1 = []
A1 = []
B2 = []
A2 = []

for d in dihedrals_rad:
    if -np.pi < d <= -1.5:
        A1.append(d) # + (np.pi-1.49))
    if -1.5 < d <= 0:
        B1.append(d) # - (np.pi-1.49))
    if 0 < d <= 1.5:
        A2.append(d) # + (np.pi-1.49))
    if 1.5 < d <= np.pi:
        B2.append(d) # - (np.pi-1.49))



fig = plt.figure(figsize=(20, 10))
ax1 = fig.add_subplot(111)
plt.xlabel('Radians')
plt.ylabel('Counts')
modes = ['A1', 'B1', 'A2', 'B2']
hist_data = []
for n, i in enumerate([A1, B1, A2, B2]):
    percent = round((len(i) / len(dihedrals_rad)) * 100, 0)
    lab = '{} = {:.2g}%'.format(modes[n], percent)
    counts, bins, _ = ax1.hist(i, bins=25, label=lab)
    for j in range(len(counts)):
        hist_data.append((modes[n], bins[j], bins[j+1], counts[j]))


# After the loop, write the data to a file
with open("histogram_data.txt", "w") as f:
    f.write("Mode,BinStart,BinEnd,Count\n")
    for mode, bin_start, bin_end, count in hist_data:
        f.write(f"{mode},{bin_start},{bin_end},{count}\n")


plt.legend(loc='upper left', fancybox=True, shadow=True, ncol=2)
plt.savefig('SimulatedTolBindingModes.png', dpi=300, bbox_inches='tight')
plt.savefig('BMAnalysisB50.pdf', bbox_inches='tight')

plt.show()

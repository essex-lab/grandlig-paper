
import pandas as pd
import os
import argparse
from openmm import *
from openmm.app import *
from pdbfixer import PDBFixer
from openff.toolkit import Molecule, Topology
from openmm.unit import *

parser = argparse.ArgumentParser()
parser.add_argument('--df')
parser.add_argument('-d', "--delim", help="Give the csv files delimiter")
parser.add_argument('-s', "--smi", help="The smiles column name")
parser.add_argument('-n', "--name", help="The id column name")


args = parser.parse_args()


# get working dir
wd = os.getcwd()
print(wd)

df = pd.read_csv(args.df, delimiter=args.delim)

for i in range(len(df)):
    x = df.iloc[i]
    id = str(x[args.name])
    print(id)

    try:
        os.mkdir(id)
    except:
        continue

    os.chdir(id)  # Move into the directory
    try:
        os.mkdir("Mu_Ex")
    except:
        continue

    os.chdir("Mu_Ex")

    for j in range(3):
        os.mkdir(f"repeat_{j+1}")

    smile = x[args.smi]
    

    lig = Molecule.from_smiles(smile)
    # Generate a conformer to be used as atomic coordinates
    lig.generate_conformers(n_conformers=1)
    lig.to_file(f"{id}.sdf", 'sdf')
    for i in lig.atoms:
        i.metadata['residue_name'] = "L02"

    top = lig.to_topology()

    top.to_file(f"{id}.pdb")

    
    # Solvate the topology
    fixer = PDBFixer(f'{id}.pdb')
    fixer.addSolvent(
        boxSize=Vec3(3, 3, 3) * nanometer, ionicStrength=0.15 * molar
    )

    with open(f"{id}_solvated.pdb", "w") as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


    os.chdir(wd)



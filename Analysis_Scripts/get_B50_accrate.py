"""
Script to find the closest simulated B value to B50 and extract the final acceptance rate. Script requires tidying up a bit.
- Will Poole
"""

# Modules to import
import argparse
import os

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from tqdm import tqdm
import glob 
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

small_font = 16
medium_font = 18
large_font = 20

plt.rc("figure", titlesize=large_font)
plt.rc("font", size=small_font)
plt.rc("axes", titlesize=small_font)
plt.rc("axes", labelsize=medium_font)
plt.rc("xtick", labelsize=small_font)
plt.rc("ytick", labelsize=small_font)
plt.rc("legend", fontsize=small_font)

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument(
    "--B50s", help="Input mu file to egt muex from", type=str, default=None
)
parser.add_argument(
    "--smis", help="path to a file mapping ligand identifiers to SMILES (columns: Ligand/Name and SMILES/Smiles)", 
    type=str, default=None
)

args = parser.parse_args()


def get_acc(file):
    with open(file, "r") as f:
        lines = f.readlines()
        final_line = lines[-1].split()
        acc_rate = final_line[9].strip('(')
    return float(acc_rate)

B50s = pd.read_csv(args.B50s, sep="\t", converters={"Ligand": str, "Name": str})
smis = pd.read_csv(args.smis, delim_whitespace=True, converters={"Ligand": str, "Name": str})
smis.rename(columns={"Name": "Ligand", "SMILES": "Smiles"}, inplace=True)
# Merge on name columsn
df = B50s.merge(smis, on="Ligand")
df = df[["Ligand", "B50", "Smiles"]]
mols = [Chem.MolFromSmiles(smi) for smi in df["Smiles"]]
hacs = [m.GetNumHeavyAtoms() for m in mols]
n_rot = [rdMolDescriptors.CalcNumRotatableBonds(m) for m in mols]

print(df)


ligands = df["Ligand"].to_list()
B50s = df["B50"].to_list()
import re

mean_acc_per_ligand = []
raw_acc_per_ligand = []
sem_acc_per_ligand = []
closest_b_per_ligand = []
for i in range(len(ligands)):
    lig_name = ligands[i]
    lig_name_glob = re.sub("([\[\]])", "[\\1]", lig_name)
    lig_done = False
    b50 = B50s[i]
    acc_rate_per_repeat = []
    for repeat in [1, 2, 3, 4]:
        closest_file, closest_b = None, None
        # base_paths = [
        #     f"/home/will/data_2/Grand/T4L99A/Titrations/{lig_name_glob}/Apo/150ps/repeat_{repeat}/*.*/T4L99A.log",
        #     f"/home/will/data_5/Grand/T4L99A/Titrations_newligands/{lig_name_glob}/Apo/150ps/repeat_{repeat}/*.*/gcncmc.log",
        #     f"/home/will/data_6/Grand/T4L99A/Titrations_newligands/{lig_name_glob}/Apo/150ps/repeat_{repeat}/*.*/gcncmc.log",
        # ]
        base_paths = [
            f"/home/will/data_3/Grand/Occluded_Tests/MUP1/Simulations/{lig_name_glob}/repeat_{repeat}/*.*/MUP1NCMC.log",
        ]

        found_files = []
        found_bs = []
        for path in base_paths:
            matched_files = glob.glob(path)  # Find all matches
            for file in matched_files:
                b = float(file.split("/")[-2])
                if closest_b is None or abs(b - b50) < abs(
                    closest_b - b50
                ):
                    closest_b = b
                    closest_file = file

        if closest_file:
            print(f"Closest file: {closest_file}")
            try:
                acc_rate_per_repeat.append(get_acc(closest_file))
            except:
                pass

        else:
            raise Exception(f"No matching files found for ligand {ligands[i]} repeat {repeat}.")
    raw_acc_per_ligand.append(acc_rate_per_repeat)
    closest_b_per_ligand.append(closest_b)
    mean_acc_per_ligand.append(np.mean(acc_rate_per_repeat))
    sem_acc_per_ligand.append(np.std(acc_rate_per_repeat) / np.sqrt(len(acc_rate_per_repeat)))

for i in range(len(ligands)):
    print(f"{ligands[i]}-{closest_b_per_ligand[i]}: {mean_acc_per_ligand[i]} +/- {sem_acc_per_ligand[i]}")


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Convert to NumPy arrays for easier sorting
hacs = np.array(hacs)
mean_acc_per_ligand = np.array(mean_acc_per_ligand)
sem_acc_per_ligand = np.array(sem_acc_per_ligand)
ligands = np.array(ligands)
raw_acc_per_ligand = np.array(raw_acc_per_ligand)

# Sort by HAC
sorted_indices = np.argsort(hacs)
hacs = hacs[sorted_indices]
mean_acc_per_ligand_hac = mean_acc_per_ligand[sorted_indices]
sem_acc_per_ligand_hac = sem_acc_per_ligand[sorted_indices]
ligands_hac = ligands[sorted_indices]
raw_acc_per_ligand_hac = raw_acc_per_ligand[sorted_indices]

# Set seaborn style
sns.set_style("whitegrid")

# Create multi-panel figure
fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"wspace": 0.3})

## Panel 1: Scatter Plot (Mean Acceptance Rate vs HAC)
axes[0].errorbar(
    hacs,
    mean_acc_per_ligand_hac,
    yerr=sem_acc_per_ligand_hac,
    fmt="o",
    markersize=8,
    color="royalblue",
    ecolor="crimson",
    capthick=2,
    capsize=5,
    alpha=0.8,
    label="Mean Acc Rate",
)

df = pd.DataFrame({'HACs': hacs, 'Mean Acc Rate': mean_acc_per_ligand_hac})
df.to_csv("MUP1_HACs_vs_mean_acc_rate.csv", index=False)

axes[0].set_xlabel("Heavy Atom Count (HAC)", fontsize=14, fontweight="bold")
axes[0].set_ylabel("Mean Acceptance Rate", fontsize=14, fontweight="bold")
axes[0].set_title(
    "Mean Acceptance Rate vs Heavy Atom Count", fontsize=16, fontweight="bold"
)
axes[0].tick_params(axis="both", labelsize=12)

## Panel 2: Bar Plot (Mean Acceptance Rate vs Ligand Name)
axes[1].bar(
    ligands_hac,
    mean_acc_per_ligand_hac,
    yerr=sem_acc_per_ligand_hac,
    capsize=5,
    color="royalblue",
    edgecolor="black",
    alpha=0.8,
)

df = pd.DataFrame({'Ligand': ligands_hac, 'Mean Acc Rate': mean_acc_per_ligand_hac})
df.to_csv("MUP1_ligands_vs_mean_acc_rate.csv", index=False)

# Plot raw data for each repeat as black dots
for i, raw_acc in enumerate(raw_acc_per_ligand_hac):
    x = np.full(len(raw_acc), i)  # x-coordinates for the dots
    axes[1].scatter(
        x,
        raw_acc,
        color="black",
        alpha=0.6,
        zorder=5,)

axes[1].set_xlabel("Ligand Name", fontsize=14, fontweight="bold")
axes[1].set_ylabel("Mean Acceptance Rate", fontsize=14, fontweight="bold")
axes[1].set_title("Mean Acceptance Rate vs Ligand", fontsize=16, fontweight="bold")
axes[1].tick_params(axis="y", labelsize=12)
axes[1].set_xticks(range(len(ligands)))  # Ensure correct x-tick positions
axes[1].set_xticklabels(
    ligands_hac, rotation=45, ha="right", fontsize=10
)  # Rotate with right alignment


# Optimize layout
plt.tight_layout()
plt.suptitle("MUP1")
plt.savefig("acceptance_rate_vs_hac_and_ligand.pdf", bbox_inches='tight')
plt.show()


# plt.clf()

# # Convert to NumPy arrays for easier sorting
# n_rot = np.array(n_rot)

# # Sort by Rot
# sorted_indices = np.argsort(n_rot)
# n_rot = n_rot[sorted_indices]
# print(n_rot)
# mean_acc_per_ligand_rot = mean_acc_per_ligand[sorted_indices]
# sem_acc_per_ligand_rot = sem_acc_per_ligand[sorted_indices]
# print(ligands)
# ligands_rot = ligands[sorted_indices]
# print(ligands_rot)

# # Set seaborn style
# sns.set_style("whitegrid")

# # Create multi-panel figure
# fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"wspace": 0.3})

# ## Panel 1: Scatter Plot (Mean Acceptance Rate vs HAC)
# axes[0].errorbar(
#     n_rot,
#     mean_acc_per_ligand_rot,
#     yerr=sem_acc_per_ligand_rot,
#     fmt="o",
#     markersize=8,
#     color="royalblue",
#     ecolor="crimson",
#     capthick=2,
#     capsize=5,
#     alpha=0.8,
#     label="Mean Acc Rate",
# )

# axes[0].set_xlabel("Number of Rotatable Bonds", fontsize=14, fontweight="bold")
# axes[0].set_ylabel("Mean Acceptance Rate", fontsize=14, fontweight="bold")
# axes[0].set_title(
#     "Mean Acceptance Rate vs Number of Rotatable Bonds", fontsize=16, fontweight="bold"
# )
# axes[0].tick_params(axis="both", labelsize=12)

# ## Panel 2: Bar Plot (Mean Acceptance Rate vs Ligand Name)
# axes[1].bar(
#     ligands_rot,
#     mean_acc_per_ligand_rot,
#     yerr=sem_acc_per_ligand_rot,
#     capsize=5,
#     color="royalblue",
#     edgecolor="black",
#     alpha=0.8,
# )

# axes[1].set_xlabel("Ligand Name", fontsize=14, fontweight="bold")
# axes[1].set_ylabel("Mean Acceptance Rate", fontsize=14, fontweight="bold")
# axes[1].set_title("Mean Acceptance Rate vs Ligand", fontsize=16, fontweight="bold")
# axes[1].tick_params(axis="y", labelsize=12)
# axes[1].set_xticks(range(len(ligands)))  # Ensure correct x-tick positions
# axes[1].set_xticklabels(
#     ligands_rot, rotation=45, ha="right", fontsize=10
# )  # Rotate with right alignment

# plt.show()

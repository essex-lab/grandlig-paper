"""
Main script to take in occupancy data from a GCNCMC simulation and plot a titration curve. This script is intened to be used for multiple ligands with an indivdiual
 titration curve for each ligand plotted. Curves are plotted on the B and concentration scale. A third curve is plotted after boootstrapping the data between repeats.
A numpy array with the Adams values, mean fitted curve parameters, error in curve parameters and the free enegry data for each repeat is written out to a .npy file for future processing and plotting.
"""
# Modules to import
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
import pandas as pd

# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument(
    "--mu_file",
    help="Input mu file to egt muex from",
    type=str,
    default=None
)
args = parser.parse_args()

try:
    mu_df = pd.read_csv(
        args.mu_file,
        converters={"Name": str},
    )
except:
    mu_df = pd.read_csv(
        args.mu_file, delim_whitespace=True,
        converters={"Name": str},
    )

ligands = mu_df["Name"].to_list()

# Combine all ligands and repeats into a single dataframe
combined_data = []
for i, lig in enumerate(ligands):
    file_name = lig + "_NOpost_procress_data.csv"

    # Read the CSV file without headers
    df = pd.read_csv(file_name, header=None, na_values=["nan", "NaN", " nan", " ", ""])

    # Set custom column names
    column_names = [f"B_{lig}"] + [f"r{j}_{lig}" for j in range(1, df.shape[1])]
    df.columns = column_names

    # Keep only B and r1–r4
    columns_to_keep = [f"B_{lig}"] + [f"r{j}_{lig}" for j in range(1, 5)]
    df = df[columns_to_keep]

    # Reset index to avoid alignment issues when concatenating
    df.reset_index(drop=True, inplace=True)

    combined_data.append(df)

# Concatenate all DataFrames column-wise
final_combined_df = pd.concat(combined_data, axis=1)


# Show result
print(final_combined_df.head())
# Save the combined DataFrame to a new CSV file
output_file = "combined_data.csv"
final_combined_df.to_csv(output_file, index=False)
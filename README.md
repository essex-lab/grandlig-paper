# Supplementary Material for "Accelerating Fragment Baed Drug Discovery using Grand Canonical Nonequilibrium Candidate Monte Carlo"

This repository contains the data and scripts necessary to reproduce the results presented in our manuscript, _"Accelerating Fragment Baed Drug Discovery using Grand Canonical Nonequilibrium Candidate Monte Carlo"_, by W. G. Poole, M. L. Samways, D. Branduardi, R. D. Taylor, M. L. Verdonk, and J. W. Essex (2024).

There are separate directories for the bulk concentrations simulations as well as the three test systems: beta-cyclodextrin, T4L99A, and MUP1 with the scripts and input files required to run the GCNCMC/MD simulations. Where the necessary scripts to reproduce the published data can be found. The required analysis can also be found in each sub-directory. More information on how to perform the simulations and analysis can be found within each directory. The _grandlig_ module required to run these simulations can be downloaded [here](https://github.com/essex-lab/grand-lig').

If you encounter any problems with these scripts, please open an issue or contact the authors.

# Overview
## Concentration
Input files to run the concentration simulations.

## HelperScripts
In this folder you will find general scripts which may be useful for running GCNCMC simulations in the future.

## Ligands
Text files containing the name, smile string and calculated excess chemical potentials for the ligands studied in this paper. A notebook is provided to recreate the figures in the SI depicting the structure of the molecules.

## Mu_Ex_Calculation
Example scripts and input files to calculate the excess chemical potential for some of the ligands studied.

## HostGuest
Scripts and input files to reproduce results in the paper regarding the host guest beta-cyclodextrin system.

## T4L99A
Scripts and input files to reproduce results in the paper regarding the T4L99A system.

## MUP1
Scripts and input files to reproduce results in the paper regarding the MUP1 system.

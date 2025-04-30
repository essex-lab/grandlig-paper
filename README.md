# Supplementary Material for "Accelerating Fragment Baed Drug Discovery using Grand Canonical Nonequilibrium Candidate Monte Carlo"

This repository contains the data and scripts necessary to reproduce the results presented in our manuscript, _"Accelerating Fragment Baed Drug Discovery using Grand Canonical Nonequilibrium Candidate Monte Carlo"_, by W. G. Poole, M. L. Samways, D. Branduardi, R. D. Taylor, M. L. Verdonk, and J. W. Essex (2024).

There are separate directories for the bulk concentrations simulations as well as the three test systems: beta-cyclodextrin, T4L99A, and MUP1 with the scripts and input files required to run the GCNCMC/MD simulations. The required analysis can also be found in each sub-directory. More information on how to perform the simulations and analysis can be found within each directory. The _grandlig_ module required to run these simulations can be downloaded [here](https://github.com/essex-lab/grand-lig'). Installation instructions are provided on the main repo.

If you encounter any problems with these scripts (no matter how small or large), please open an issue and it will get fixed ASAP.

# Overview
## Summary
The simplest approach to reproduce the results of the paper is to use the scripts in `GCNCMC_Simulation_Scripts` with different arguments to recreate the results. This will require the user to set up their own folder structures. A lot of simulations have been performed for this paper and therefore regenerating these results is no small task and will require the use of HPCs. 

While best efforts have been taken to ensure everything in this repo runs and works as expected, there are clearly many scripts and input files provided. If you run into any issues please open an issue and we will amend them ASAP. 

## Manifest
### Concentration
Input files to run the concentration simulations.

### Ligands
Text files containing the name, smile string and calculated excess chemical potentials for the ligands studied in this paper. A notebook is provided to recreate the figures in the SI depicting the structure of the molecules.

### Mu_Ex_Calculation
Example scripts and input files to calculate the excess chemical potential for some of the ligands studied.

### HostGuest
Scripts and input files to reproduce results in the paper regarding the host guest beta-cyclodextrin system.

### T4L99A
Scripts and input files to reproduce results in the paper regarding the T4L99A system.

### MUP1
Scripts and input files to reproduce results in the paper regarding the MUP1 system.

### Basic_Simulation
Basic scripts which can be used to run sphere based GCNCMC simulations. These scripts use `argparse` to allow the used to provide all the relevent information to run a simulation. While input files and run scripts are provided in the above sections, these scripts can also be used to reproduce the results of the paper. This route provides more flexibiliy may be easier to work with.

### GCNCMC_Simulation_Scripts
Set of scripts required to perform various pieces of analysis. Again, these can be manipulated to used as inspiration for your own analysis and beyond.

### Outputs
Collated output data which matches the publication.

# Generalisation
This repo has many many input files. This is so that the results in the paper are reproducible. To make things simpler and clearer this is a more general overview to the GCNCMC protocols.

For any given ligand and system:
1. An excess chemical potential needs to be calculated for the ligand. This is basic hydration/solvation free energy calculation. Note that as mentioned in the paper, for suffciently dilute concentrations and hydration free energy will suffice. For higher concentrations, depending on application, it may be more precise to perform a solvation free energy calculation in a solution of that concentration. Examples are found in `T4L99A/FEs/Solv`, [the main repo](https://github.com/essex-lab/grand-lig/examples), and in `Mu_Ex_Calculations`. 

2. Now the user must make choices on what kind of simulation they want to run. To find occluded sites you may choose to define a sphere which covers the whole protein, and run at a certain concentration. I.e. `T4L99A/BindingSite`. Or to calculate free energies via titration you may want to set the sphere to cover the binding site and run simulations at a range of cocnentraions/B values. I.e. `T4L99A/Titrations`
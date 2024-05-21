# Supplementary Material for "Small Molecule Binding Sites, Poses, And Free Energies with Grand Canonical Nonequilibrium Candidate Monte Carlo"

This repository contains the data and scripts necessary to reproduce the results presented in our manuscript, _"Small Molecule Binding Sites, Poses, And Free Energies with Grand Canonical Nonequilibrium Candidate Monte Carlo"_, by W. G. Poole, M. L. Samways, D. Branduardi, M. L. Verdonk, and J. W. Essex (2022).

There are separate directories for the bulk concentrations simulations as well as the two test systems: beta-cyclodextrin, and T4L99A with the scripts and input files required to run the GCNCMC/MD simulations. For T4L99A, the directory is further sub-divided into ```Binding Site/```, ```Binding Modes/``` and Titrations, where the necessary scripts to reproduce the published data can be found. The required analysis can also be found in each sub-directory. More information on how to perform the simulations and analysis can be found within each directory. The _grand_ module required to run these simulations can be downloaded [here](https:github.com/essex-lab/ligand-grand/').

If you encounter any problems with these scripts, please contact the authors.

# Running these scripts:

## T4L99A Binding site:

Note: The script is set to run for 50 moves. This is less than the paper but should be suffcient to find the binding site and run in appox. 1 hour on a gtx1080.

```
cd T4L99A/BindingSite
mkdir repeat_1
cd repeat_1
python ../Benzene_Binding_site.py -p ../T4L99A_equiled.pdb -l ../Benzene.pdb -x ../Benzene.xml -r L02 -c 0.5
```

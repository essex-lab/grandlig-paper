# Folder for MuEx calculations

A pre-requisite for GCNCMC calculations is the excess chemical potential (mu_ex) of the molecule being simualated. This can be calculated by performing a basic hydration free energy calculation of your molecule. 

In this folder we have provided a quick way of setting up and calculating these for all the ligands studied in work.

To generate coordinates and folder structures for e.g. T4:

```
python part_1_setup_openff.py --df t4_ligands.csv -d ';' -s 'Smiles' -n "Name"
```
This will make a folder structure like: `Ligand/Mu_Ex/repeat_?`

To run a simulation e.g.:
```
cd Benzene/Mu_Ex/repeat_1
python ../../../calcMuEx.py -p ../Benzene_solvated.pdb -l ../Benzene.sdf
```

In these particular examples, the mu ex is calculated in one single script which loops through a certain number of lambda windows in a single simulation. 

To make use of parallelisation, each individual lambda could be run as a seperate simulation with the results combined at the end to calculate a final free energy. An example of this is can be found in `T4L99A/FEs/Solv`.




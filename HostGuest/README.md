# Host Guest GCNCMC Simulations.

In this folder, you should find everything rqeuired to run host guest GCNCMC simulations.

## Files and Folders
- `Ligand_pdbs_xmls/`: Folder containing all ligand 3D coordinates (.pdb) and GAFF/AM1BCC forcefield parameters (.xml)
- `Titrations/`: A starting place to run GCNCMC Host Guest titrations.
- `Titrations/HG_B_values.csv`: A file containing the simulated B values for each ligand
- `bCD_FINAL.xml`: Parameters for the beta-cyclodextrin host
- `bCD_Equiled.pdb`: Equilibrated coordinates for the bCD host

## GCNCMC
To run a single simulation for a particular ligand at a given B value you can run:

```
python basic_NCMC_sim_B_HG.py --pdb bCD-Equiled.pdb  --lig LIG_PDB --lxml LIG_XML --B BVALUE --st 50 --sph_atoms C34 C7 --rad 5 --nmoves 2000 --mdpermove 500
```

Where `LIG_PDB` and `LIG_XML` files are in `Ligand_pdbs_xmls` and `BVALUE` can be found in `Titrations/HG_B_values.csv`. 

As an example using `benzaldehyde` at a B value of `-13.24`:

```
mkdir -p Titrations/benzaldehyde/repeat_1/-13.24/
cd Titrations/benzaldehyde/repeat_1/-13.24/

python ../../../../../GCNCMC_Simulation_Scripts/basic_NCMC_sim_B_HG.py --pdb ../../../../bCD-Equiled.pdb --hxml ../../../../bCD_FINAL.xml --lig ../../../../Ligand_pdbs_xmls/benzaldehyde.pdb --lxml ../../../../Ligand_pdbs_xmls/benzaldehyde.xml --B -13.24 --st 50 --sph_atoms C34 C7 --rad 5 --nmoves 2000 --mdpermove 500
```

To run a full titration calculation you will need to perform simulations at multiple B values.

For analysis, see the `Outputs/HostGuest` folder.

## FEP
Privided are representative frames from GCNCMC simulations of both the primary and secondary poses of each HG pair - these are the starting structures for FEP simulations. 

To perform a FE calculation (example):

```
cd FEP/1-methylcyclohexanol/Complex_Prim
mkdir repeat_1; cd repeat_1

for i in $(seq 0 40); do mkdir $i; cd $i; python HostGuestComplex_free_energy.py -p ../../Primary.pdb -x ../../*.xml -r L02 848 -lam $i; cd ..; done

```

Each lambda value is run indivdually with a numpy array written out - `U_matrix_i.npy` of shape `(40, 40, 1000)`. This data can be combined to give the full matrix for each repeat:

```
cd FEP/1-methylcyclohexanol/Complex_Prim/repeat_1
python Free_Energy_Scripts/Combine_FE_numpy_arrays.py -u */*.npy
```

Analysis can then be performed on this combined matrix. See `Outputs/HostGuest/FEP` for examples.
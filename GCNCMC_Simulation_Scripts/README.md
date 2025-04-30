# Basic Simulation Scripts

These scripts are provided here for convinience and can also be found in the main repo (https://github.com/essex-lab/grand-lig/tree/master/examples/Basic_Simulation). 

They can be used, with some manipulations, to reproduce the results of this paper as follows:

## Manifest
- `basic_NCMC_sim_B_HG.py` a simple script to perform a GCNCMC simulation at a user defined B value for the Host Guest system. 
- `basic_NCMC_sim_B.py` a simple script to perform a GCNCMC simulation at a user defined B value for protein systems.
- `basic_NCMC_sim_mu.py` a simple script to perform a GCNCMC simulation at a user defined excess chemical potential and concentration for protein systems - the B value is calculated under the hood.


## Examples
### Host Guest
### Titrations
For any given ligand and B value, a simulation can be run as:
```
python basic_NCMC_sim_B_HG.py --pdb bCD-Equiled.pdb  --lig LIG_PDB --lxml LIG_XML --B BVALUE --st 50 --sph_atoms C34 C7 --rad 5 --nmoves 2000 --mdpermove 500
```
Input files are found in `../HostGuest` and B values for each ligand are provided in a csv file in `../HostGuest/Titrations/HG_B_values.csv` 

### T4L99A
### Site Finding using 0.5 M Benzene
Input files and a script for these simulations are found at `../T4L99A/BindingSite`

Alternatively:
```
python basic_NCMC_sim_mu.py --pdb T4L99A_equiled.pdb --lig Benzene.pdb --xml Benzene.xml --mu -0.679 --conc 0.5 --sph_resns PHE GLU --sph_resis 105 12 --rad 26.5 --hmr --nmoves 100 --mdpermove 4000
```

### Toluene Binding Modes
Note: This simulation may require a lot of moves. Input files: `../T4L99A/BindingModes`

```
python basic_NCMC_sim_B.py --pdb T4NoCoSolv.pdb --lig Methylbenzene.pdb --xml Methylbenzene.xml --B -7.3408 --sph_resns LEU ALA --sph_resis 85 100 --rad 8 --hmr --nmoves 5000 --mdpermove 500
```

### T4 Titrations
For any given ligand and B value, a simulation can be run as:
```
python basic_NCMC_sim_B.py ---pdb T4NoCoSolv.pdb --lig LIG_PDB --xml LIG_XML --B BVALUE --sph_resns LEU ALA --sph_resis 85 100 --rad 8 --hmr --nmoves 1500 --mdpermove 500
```

Input files are found in `../T4L99A/Titrations` and B values for each ligand are provided in a csv file.

### MUP1
### Site Finding using 0.5 M of ligand 07
Input files and a script for these simulations are found at `../MUP1/BindingSite`

```
python basic_NCMC_sim_mu.py --pdb EquiledMUP1_with07.pdb --lig 07.pdb --xml 07.xml --mu -2.51 --conc 0.5 --sph_resns GLY --sph_resis 119 --rad 22 --hmr --nmoves 1000 --mdpermove 5000
```

### MUP1 Titrations
For any given ligand and B value, a simulation can be run as:
```
python basic_NCMC_sim_B.py ---pdb ProteinEquiled.pdb --lig LIG_PDB --xml LIG_XML --B BVALUE sph_resns PHE LEU --sph_resis 57 106 --rad 5.5 --hmr --nmoves 2000 --mdpermove 1500
```

Input files are found in `../MUP1/Titrations` and B values for each ligand are provided in a csv file.




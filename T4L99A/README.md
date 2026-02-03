
# T4L99A GCNCMC Simulations

In this folder, you should find everything rqeuired to run host guest GCNCMC simulations.

## Files and Folders
- `Ligand_pdbs_xmls/`: Folder containing all ligand 3D coordinates (.pdb) and GAFF/AM1BCC forcefield parameters (.xml)
- `BindingSite/`: A starting place for mixed solvent MD simulations
- `BindingMode/`: A starting place for investigating Toluene/T4 binding modes

- `Titrations/`: A starting place to run GCNCMC T4L99A titrations.
- `Titrations/T4_B_values.csv`: A file containing the simulated B values for each ligand

# Running these scripts:

### T4L99A Binding site:

Note: The script is set to run for 50 moves. This is less than the paper but should be suffcient to find the binding site and run in appox. 1 hour on a gtx1080.

```
cd T4L99A/BindingSite
mkdir repeat_1
cd repeat_1
python ../Benzene_Binding_site.py -p ../T4L99A_equiled.pdb -l ../Benzene.pdb -x ../Benzene.xml -r L02 -c 0.5
```

In short, after a few moves, one would expect to see a Benzene molecule bound in the small, hydrophobic binding pocket. Seen by aligning the xtal structure, 181l, to the simulated trajecotry in pymol.

```
pymol T4L99AWithghosts.pdb

# In pymol
load T4L99A.dcd
fetch 181l  # Requires internet connection
align 181l T4L99AWithghosts.pdb
```

To analyse:

The following combined the final 5ns from each repeat, aligns the trajectory and perfoms a grid analysis

```
python cat_last_5ns.py -t repeat_?/gcncmc.dcd -p  repeat_1/Protein_Ghosts.pdb

python aligntraj.py -di combined-new.dcd -p  repeat_1/Protein_Ghosts.pdb -do aligned.dcd -r ../../aligntoxtal.pdb

python grid_analysis_all_atom_basicCosolv_cuda.py -di aligned.dcd -p repeat_1/Protein_Ghosts.pdb -r L02 -dt 10 -t 0 -o Benzene_grid.cube
done
```

### T4L99A Binding mode:
```
cd T4L99A/BindingModes
mkdir repeat_1
cd repeat_1
python ../Toulene_BM_simulations.py
```

### Titration Calculations
As an example using `Benzofuran` at a B value of `-9.58`:

```
mkdir -p Titrations/Benzofuran/repeat_1/-9.58/
cd Titrations/Benzofuran/repeat_1/-9.58/

python ../../../../../GCNCMC_Simulation_Scripts/basic_NCMC_sim_B.py --pdb ../../../../Native_APO.pdb --lig ../../../../Ligand_pdbs_xmls/Benzofuran_H.pdb --xml ../../../../Ligand_pdbs_xmls/Benzofuran_H.xml --B -9.58 --st 150 --sph_resns LEU ALA --sph_resis 85 100 --rad 8 --hmr --nmoves 1500 --mdpermove 500
```

To run a full titration calculation you will need to perform simulations at multiple B values.

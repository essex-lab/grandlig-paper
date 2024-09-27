
# Running these scripts:

## MUP1 Binding site / Enhanced Cosolvent:

GCNCMC:
```
cd BindingSite
cd GCNCMC
python GCNCMCMD_forCoSolv.py --pdb ../EquiledMUP1_with07.pdb --lig ../07.pdb --xml ../07.xml --mu -2.51 --conc 0.5 --st 50 --sph CA GLY 119 --rad 22 --hmr
```

MD:
```
cd BindingSite
cd MD
basic_MD_forCoSolv.py -p ../EquiledMUP1_with07.pdb -x ../07.xml
```
To analyse:

The following combined the final 5ns from each repeat, aligns the trajectory and perfoms a grid analysis

```
python cat_last_5ns.py -t repeat_?/gcncmc.dcd -p  repeat_1/Protein_Ghosts.pdb

python aligntraj.py -di combined-new.dcd -p  repeat_1/Protein_Ghosts.pdb -do aligned.dcd -r ../../aligntoxtal.pdb

python grid_analysis_all_atom_basicCosolv_cuda.py -di aligned.dcd -p repeat_1/Protein_Ghosts.pdb -r L02 -dt 10 -t 0 -o 07_grid.cube
done
```

## Titration Calculations
The folder structure is already setup with a python file to run in each folder. To run a titration for e.g. ligand 01

```
cd Titration/01/repeat_1
```

You will then need to `cd` into each directory corresponding to an individual B value and run the simulation. You will need a full set of simulations at multiple B values to plot a titration curve. 



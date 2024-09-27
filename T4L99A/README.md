
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
cd T4L99A/BindingSite
mkdir repeat_1
cd repeat_1
python ../Benzene_Binding_site.py -p ../T4L99A_equiled.pdb -l ../Benzene.pdb -x ../Benzene.xml -r L02 -c 0.5
```

## Titration Calculations
The folder structure is already setup with a python file to run in each folder. To run a titration for e.g. Benzene

```
cd Titration/Benzene/repeat_1
```

You will then need to `cd` into each directory corresponding to an individual B value and run the simulation. You will need a full set of simulations at multiple B values to plot a titration curve. 



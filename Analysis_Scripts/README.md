# Analysis Scripts

Here are various analysis scripts required to analyse the GCNCMC simulations. Note these are minimal scripts and you may prefer to build your own. 

- `aligntraj.py` Simple script to align a simulation trajectory to the first frame, or optioanlly a reference PDB. Arguments can be found using the `-h` flag.

- `trajcat.py` Concatenates a list of trajectories into one larger trajecotry file. Useful for combining data from multiple repeats such as in MSMD simulations.

- `functions.py` Contains all useful functions pertaining to titration analysis

- `p-xylene_val111_dihedral.py` Script which reads in a simualtion trajectory and plots the Val111 dihedral angle.



- `grid_analysis_all_atom_basicCosolv_cuda.py` - Basic script to run grid analysis. Other scripts exist and the same analysis can be perfomed in VMD. 

- 

## Titrations
For titration calculations, to write out the final average occupancies of simulations of a particular ligand e.g. host-guest ligand trans-4-methylcyclohexanol:

```
python Write_avg_occs.py -f repeat_? -fn betaCD.log -frag trans-4-methylcyclohexanol
```

To plot titration curves and return a free energy e.g:
```
python Plot_Titration_curve_Multiple_lig_V2.py -G HG_ligands_mu.txt -r 5 -pro "beta CD"
```

Where `HG_ligands_mu.txt` is a comma seperated file containing the excess chemical potentials of the ligands.
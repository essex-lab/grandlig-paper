# Analysis Scripts

Here are various analysis scripts required to analyse the GCNCMC simulations. Note these are minimal scripts and you may prefer to build your own. The figures in the paper with which these scripts align with are indicated.

- `aligntraj.py` Simple script to align a simulation trajectory to the first frame, or optioanlly a reference PDB. Arguments can be found using the `-h` flag.

- `functions.py` Contains all useful functions pertaining to titration analysis

- `grid_analysis_all_atom_basicCosolv_cuda.py` - Basic script to run grid analysis. Other scripts exist and the same analysis can be perfomed with a GUI in VMD. (Figures 3 and 4)

- `p-xylene_val111_dihedral.py` Script which reads in a simulation trajectory and plots the Val111 dihedral angle. (Figure S22)

- `Plot_Combined_TitCurve_with_bar.py` Script which takes in curve parameters for each ligand and plots the titration curves on a single plot alongside a bar chart containing in the final free energy estimates. (Figures 8, S12 and S24)

- `plot_tit_curves_with_convergance.py` Script which reads in occupancy data as a function of number of moves performed and calculates the free energy estimate as over time. A convergance plot is returned. (Figure S17)

- `Plot_Titration_curves_B_C_and_Boot.py` Reads in the final occupancy data for a ligand and plots the titration curve on the B and log(conc) scale. The data is written out for further plotting and processing. 

- `trajcat.py` Concatenates a list of trajectories into one larger trajecotry file. Useful for combining data from multiple repeats such as in MSMD simulations.

- `Write_avg_occs_function_nmoves.py` Reads in the log files for a GCNCMC simulation and writes out the average occupancy across repeats as a function of the number of moves performed. 

- `Write_avg_occs.py` Reads in log files from GCNCMC simulations and writes out the final average occupancy. 

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
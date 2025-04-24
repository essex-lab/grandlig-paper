# T4L99A Outputs

- `BindingSite` : Contains only the final frame from one repeat and the `.cube` file from grid analysis. Trajectories avaliable upon request as they are too large for this repo.

- `Titrations`: Contains titration bits

- `Titrations/raw_log_files`: Contains log files from the outputs of all titration calculations.

- `Titrations/Final_average_occs`: Contains final average occupancies from all titration calculations.

- `Titrations/Convergance_plots`: Contains convergance plots for all titration calculations

- `Titrations/Titration_ResultsNEW.txt`: Final free energy estimates for each titrated ligand




## Titration Analysis
A guide to analysing these titration calculations:

1) Extract final average occupancies from the log file (Final_average_occs):

```
cd Titrations/raw_log_files
for l in *; do cd $l/Apo/150ps; python Write_avg_occs.py -f repeat_? -fn gcncmc.log -frag $l; cd ../../../; done
```

2) Calculate free energy estimate from final occupancies for all ligands:
```
python ../../../Analysis_Scripts/Plot_Titration_curves_B_C_and_Boot.py -r 8 -n 1000 --mu_file ../Mu_ex.txt
```

For convergance analysis:
```
cd Titrations/raw_log_files
for l in *; do python ../../../../Analysis_Scripts/Write_avg_occs_function_nmoves.py -f $l/Apo/150ps/repeat_? -fn gcncmc.log -frag $l; done
mv *.pkl ../Convergance_plots/
cd ../Convergance_plots
python ../../../../Analysis_Scripts/plot_tit_curves_with_convergance.py --mu_file ../../Mu_ex.txt
```


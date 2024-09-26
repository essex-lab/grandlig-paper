# Concentration Simulations

In this folder are the scripts to run concentration simulations in a simple box. Water and ligand moves are performed to acheive a desired concentration of ligand.

## Acetone
```
cd Acetone
cd From_0M
python run_bulk_conc.py
```
To exted the runs:
```
python cont_run.py -n 2
```

## Pyrimidine

Same as above.

## Analysis
The output from these simulations that we need is the log files. In these log files, the number of 'on' molecules in the system is reported at each frame (after every attempted move).  

To write out the concentration at each frame parse the gcncmc log files to `write_concs.py`
```
python write_concs.py -d repeat_?/bulk_conc_0.1.log repeat_?/bulk_conc_0.1-2.log -n 0 -v 62.39127066552801 -o Acetone_From0M.csv -r 8
```

where -d is a list of all the log files, -n is the initial number of molecules in the system 0 if starting from 0M, -v is the system volume and can be taken from the openmm simulation log -r is the number of repeats.


To plot:

```
python ../../Concentration/plot_concs_new.py -i Acetone_From0M.csv -c 0.5 -e 5000 --title "Acetone from 0.5M" --yzero -o 0.5_acetone_from_0
```

The published results are found in `../Outputs/Concentration`
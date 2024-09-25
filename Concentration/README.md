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
To write out the concentration at each frame parse the gcncmc log files to `write_concs.py`

To plot `plot_concs_new.py`
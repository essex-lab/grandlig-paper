# Concentration Simulations.

These files contain the raw data needed to reproduce the plots in Figure 2.

They are the the output of the script `../../Concentration/write_concs.py` which extracts the number of molecules in the system from a GCNCMC log file and converts to a concentration.

The plots can be produced by e.g.:

```
python ../../Concentration/plot_concs_new.py -i Acetone_From0M.csv -c 0.5 -e 5000 --title "Acetone from 0.5M" --yzero  -o 0.5_acetone_from_0
```
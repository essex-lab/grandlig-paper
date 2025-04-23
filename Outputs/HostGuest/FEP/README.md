# Host Guest FEP

Provided are combined numpy arrays of the collected potential energy sample for each `ligand/bindinmode/repeat`

To calculate the free energies for each:
```
cd Raw
for l in *; do cd $l/Complex_Prim; python ../../../../../../Free_Energy_Scripts/HostGuest/calc_eq_dG.py >> ../../../Processed/primary_fes.txt; cd ../../; done

for l in *; do cd $l/Complex_Sec; python ../../../../../../Free_Energy_Scripts/HostGuest/calc_eq_dG.py >> ../../../Processed/secondary_fes.txt; cd ../../; done

```

To combine the binding modes:

```
cd ../Processed
python ../../../../Free_Energy_Scripts/HostGuest/Calc_Total_FEP.py
```


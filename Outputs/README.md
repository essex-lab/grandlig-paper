# Titrations

To plot titration curves:

```
python ../../../Analysis_Scripts/Plot_Titration_curve_Multiple_lig_V2.py -G HG_ligands_mu.txt -r 5 -pro "beta CD"
```

This script will also output a `.pkl` file for easier plotting the second time:

```
python ../../../Analysis_Scripts/Plot_Curves_From_File.py -pkl titration_plot_data.pkl -pro beta_cd -o TEST
```


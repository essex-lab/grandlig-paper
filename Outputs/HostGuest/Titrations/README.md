# Host Guest Titrations 

## Reproducing the plots

To extract free energy estimates and reproduce the figures in the paper:

```
cd Raw

for l in *; do python -i ../../../../Analysis_Scripts/Write_avg_occs.py -f $l/repeat_? -fn betaCD.log -frag $l; done

mv *NO*.csv ../Processed

cd ../Processed

python ../../../../Analysis_Scripts/Plot_Titration_curves_B_C_and_Boot.py -r 5 -n 500 --mu_file ../HG_ligands_mu.txt

python ../../../../Analysis_Scripts/Plot_Combined_TitCurve_with_bar.py -mu ../HG_ligands_mu.txt -fe Titration_ResultsNEW.txt -rad 5 -pro 'beta-cyclodextrin' -out betaCD_titrations_withBAR

cd ../

python ../../../../Analysis_Scripts/yx_plots/PlotxvsyFreeEnergy_new.py -x AllExp.txt -y Processed/Titration_ResultsNEW.txt -x_label 'Experimental' -y_label 'GCNCMC Titration' -o HG_exp_tit

python ../../../../Analysis_Scripts/yx_plots/PlotxvsyFreeEnergy_new.py -x flat_bottom_fep_JACOBIAN.txt -y Processed/Titration_ResultsNEW.txt -x_label 'FEP' -y_label 'GCNCMC Titration' -lims -6.55 1.87 -o HG_fep_tit

python ../../../../Analysis_Scripts/yx_plots/PlotxvsyFreeEnergy_new.py -x AllExp.txt -y flat_bottom_fep_JACOBIAN.txt -x_label 'Experimental' -y_label 'FEP' -lims -6.55 1.87 -o HG_exp_fep
```
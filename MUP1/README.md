
# Running these scripts:

### MUP1 Binding site / Enhanced Cosolvent:

GCNCMC:
```
cd BindingSite
cd GCNCMC
python GCNCMCMD_forCoSolv.py --pdb ../EquiledMUP1_with07.pdb --lig ../07.pdb --xml ../07.xml --mu -2.51 --conc 0.5 --st 50 --sph CA GLY 119 --rad 22 --hmr
```

MD:
```
cd BindingSite
cd MD
basic_MD_forCoSolv.py -p ../EquiledMUP1_with07.pdb -x ../07.xml
```
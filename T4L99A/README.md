
# Running these scripts:

### T4L99A Binding site:

Note: The script is set to run for 50 moves. This is less than the paper but should be suffcient to find the binding site and run in appox. 1 hour on a gtx1080.

```
cd BindingSite
cd GCNCMC
python Benzene_Binding_site.py -p ../T4L99A_equiled.pdb -l ../Benzene.pdb -x ../Benzene.xml -r L02 -c 0.5
```

In short, after a few moves, one would expect to see a Benzene molecule bound in the small, hydrophobic binding pocket. Seen by aligning the xtal structure, 181l, to the simulated trajecotry in pymol.

```
pymol T4L99AWithghosts.pdb

# In pymol
load T4L99A.dcd
fetch 181l  # Requires internet connection
align 181l T4L99AWithghosts.pdb
```

### T4L99A Binding mode:
```
cd BindingModes
mkdir repeat_1
cd repeat_1
python ../other_one.py
```

### T4L99A Titrations



# Running these scripts:

### T4L99A Binding site:

Note: The script is set to run for 50 moves. This is less than the paper but should be suffcient to find the binding site and run in appox. 1 hour on a gtx1080.

```
cd T4L99A/BindingSite
mkdir repeat_1
cd repeat_1
python ../Benzene_Binding_site.py -p ../T4L99A_equiled.pdb -l ../Benzene.pdb -x ../Benzene.xml -r L02 -c 0.5
```

In short, after a few moves, one would expect to see a Benzene molecule bound in the small, hydrophobic binding pocket. Seen by aligning the xtal structure, 181l, to the simulated trajecotry in pymol.

```
pymol T4L99AWithghosts.pdb

# In pymol
load T4L99A.dcd
fetch 181l  # Requires internet connection
align 181l T4L99AWithghosts.pdb
```

### T4L99A Binding mode:
```
cd T4L99A/BindingSite
mkdir repeat_1
cd repeat_1
python ../Benzene_Binding_site.py -p ../T4L99A_equiled.pdb -l ../Benzene.pdb -x ../Benzene.xml -r L02 -c 0.5
```


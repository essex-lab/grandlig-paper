with open("rest_params.txt", "r") as fi:
    for line in fi.readlines():
        lig, args = line.split(":")
        cmd = f"/scratch/wp1g16/miniconda3/envs/OMM8.0/bin/python $SCRIPT_FOLDER/solvate_equilibrate.py --pdb ../../GettingRestraintParameters/{lig}/ClosestRestraintFrame.pdb --xml ../../GettingRestraintParameters/{lig}/*.xml {args}"
        print(lig, args)
        print(cmd)
        with open("run_equil_iridis.sh", "r") as ft:
            with open(f"{lig}/run_equil_iridis.sh", "w") as fo:
                for l in ft.readlines():
                    fo.write(l)
                    if l.startswith("SCRIPT_FOLDER="):
                        fo.write(f"\n{cmd}")

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to align the protein
Accepts the system PDB file and the DCD simulation data and returns a new dcd file for the aligned protein
- Will Poole
"""


import argparse
import MDAnalysis, MDAnalysis.analysis.align
import numpy as np 
import copy
# Arguments to be input
parser = argparse.ArgumentParser()
parser.add_argument("-di", "--dcdin", help="Input .dcd Simulation File", default='system.dcd')
parser.add_argument("-p", "--pdb", help="Input PDB file", default='system.pdb')
parser.add_argument("-do", "--dcdout", help="Output DCD file for the aligned protein", default='alignedsystem.dcd')
parser.add_argument("-r", "--ref", help="If you would like, suppy a pdb file to use as the reference structure to "
                                        "align the protein.", default=None)
args = parser.parse_args()

if args.pdb.split(".")[-1] == 'cif':
    from simtk.openmm.app import pdbxfile
    structure = pdbxfile.PDBxFile(args.pdb)
    u = MDAnalysis.Universe(structure, args.dcdin)
else:
    u = MDAnalysis.Universe(args.pdb, args.dcdin)  # Creates the MDAnalysis universe for the simulation data


if args.ref:
    if args.ref.split(".")[-1] == 'cif':
        structure = pdbxfile.PDBxFile(args.ref)
        u_ref = MDAnalysis.Universe(structure)
    else:
        u_ref = MDAnalysis.Universe(args.ref)
else:
    if args.pdb.split(".")[-1] == 'cif':
        from simtk.openmm.app import pdbxfile
        structure = pdbxfile.PDBxFile(args.pdb)
        u_ref = MDAnalysis.Universe(structure)
    else:
        u_ref = MDAnalysis.Universe(args.pdb)


#protein = u_ref.select_atoms('protein')
#protein.unwrap(compound='fragments')
#
#protein_center = protein.center_of_mass(wrap=True)
#dim = u_ref.dimensions
#box_center = np.sum(dim, axis=0) / 2
#u_ref.atoms.translate(box_center - protein_center)
#
#not_protein = u_ref.select_atoms('not protein')
#
#not_protein.wrap(compound='residues')

# Details te alignment parameters
a = MDAnalysis.analysis.align.AlignTraj(u, u_ref, select="name CA", filename=args.dcdout)
a.run()  # Runs the alignment

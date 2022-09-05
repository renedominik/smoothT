A software for the construction of simple transition pathways from PDB ensembles. 
Inputs:
- a list of pdb files
- a starting conformation (PDB)
- a final conformation (PDB)
- an energy identifier (each PDB file contains a line with this identifier-string followd by the energy value - space separated)
- a RMSD cutoff, defining conformation to be neighbors for graph construction


COMPILATION:

 g++ -O3 smoothT.cpp -o smoothT

TUTORIAL:
 
 smoothT -h


BASIC EXAMPLE:

 smoothT -l LIST_OF_PDBS.txt -b START.pdb -e FINAL.pdb -i ENERGY_IDENTIFIER -o OUTDIR -d RMSD_CUTOFF 
 
 
ADVANCED FEATURES:

- select atoms used for RMSD calculation by chains and types 
- smoothT can use either a list of identical molecules, or one can provide alignments to define atom positions to be used

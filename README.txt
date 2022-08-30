A software for the construction of simple transition pathways.
Input are:
- a list of pdb files
- a starting conformation
- a final conformation
- an energy identifier (each PDB file contains a line with this identifier and the energy value)
- a RMSD cutoff, defining conformation to be neighbors

COMPILATION:

 g++ -O3 smoothT.cpp -o smoothT

TUTORIAL:
 
 smoothT -h

BASIC EXAMPLE:
 smoothT -l LIST_OF_PDBS.txt -b START.pdb -e FINAL.pdb -i ENERGY_IDENTIFIER -o OUTDIR -d RMSD_CUTOFF 
 
 

# pdb2S2
pdb2S2 computes S2 order parameters from molecular dynamics trajectories or conformational ensembles 

pdb2S2 program takes in input 4 arguments:

- the file name of a molecular dynamics trajectory in PDB format (each snapshot separated by MODEL/ENDMDL records).

- the name of the first atom of the vector for which S2 is computed. If this is located on the preceding residue the name is preceded by a "-".

- the name of the second atom of the vector for which S2 is computed. If this is located on the following residue the name is preceded by a "+".

- the name of the output file

The output file will list the S2 order parameter for each residue containing the vector.

S2 is computed according to the equation (using second order spherical harmonics) reported in the book:

Protein NMR Spectroscopy
Principles and Practice
Authors: John Cavanagh, Nicholas J. Skelton, Wayne J. Fairbrother, Mark Rance, Arthur G. Palmer III
Academic Press

CONTACT:  

Federico Fogolari  
Dipartimento di Scienze Matematiche, Informatiche e Fisiche  
Universita' di Udine  
Via delle Scienze 206  
33100 Udine - Italy  
Tel ++39 0432 494320    
E-mail federico.fogolari@uniud.it  

COMPILATION:

The program is compiled with: 

cc pdb2S2.c -o pdb2S2 -lm

RUNNING pdb2S2

Usage examples using the file sample.pdb provided here.

Computing S2 for the vector N-HN:
./pdb2S2 sample.pdb N H sample.S2 

Computing S2 for the vector C-N where C is on the preceding residue:
./pdb2S2 sample.pdb -C N sample.S2 

Computing S2 for the vector C-N where N is on the following residue:
./pdb2S2 sample.pdb C +N sample.S2 

OUTPUT

sample.S2 (e.g. for the command ./pdb2S2 sample.pdb N HN sample.S2) will contain the lines:

GLN 2 A 0.406300
ARG 3 A 0.854538
THR 4 A 0.801300
LYS 6 A 0.781181
ILE 7 A 0.820819
GLN 8 A 0.887599
VAL 9 A 0.828375
TYR 10 A 0.886492
SER 11 A 0.881635
......
......
......


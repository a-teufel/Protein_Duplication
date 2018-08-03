#script to translate pdb files int csv alignments, orginally written by M. Johnson, updated by A. Teufel.
#Some analysis scripts require the data to be in this format

#!/anaconda/bin/python

import glob
from Bio.PDB import *


def extract_sequences(file):
    p=PDBParser()
    structure=p.get_structure('X', file)
    #print(file, end = " ")
    of.write(file + ", ")
    for model in structure:
        for chain in model:
            ppb=PPBuilder()
            for pp in ppb.build_peptides(chain):
                of.write(str(pp.get_sequence()) + ", " )
    of.write("\n")


for i in range(1, 500, 10):
    
	of = open("M_SIM1_%s.csv" %i, "w")

	for file in glob.glob("rep*/%s.pdb" %i):
    		extract_sequences(file)
    
	of.close()

# to run this code from command line: ./extract_sequence.py


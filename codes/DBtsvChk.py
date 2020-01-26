"""
   ::Purpose::
   Check if all PDBs are download and write the list of not download
   This code is optimized for O2 Slurm cluster
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 13, 2018"


import jctools
import pandas as pd
from Bio import SeqIO
from Bio import PDB
import os
import argparse

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help='A file list PDB ID that can be searched', action='store', required=True)
parser.add_argument("-o", "--output", help='Output path to store the results', action='store', required=True)
parser.add_argument("-t", "--target", help='Target directory to be searched', action='store', required=True)

args = parser.parse_args()

pdbFile = args.input
outFile = args.output
tgPath = args.target

#------------------------------------------------------------------------------
# Search PDB and extract corresponding Chains
#------------------------------------------------------------------------------
dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})
    

PDBIDs = []
for pdb_id in dfPDB['PDB']:
    fname = '{}/{}.tsv'.format(tgPath, pdb_id)
    if not os.path.isfile(fname):
        PDBIDs.append(pdb_id)

outData = pd.DataFrame(PDBIDs, columns=['PDB'])
outData.to_csv(outFile, index=False)

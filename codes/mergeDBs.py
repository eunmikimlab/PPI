"""
   ::Purpose::
   Merging PDB and chains extracted from multiple DBs
   and writing them into a file
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 16, 2018"


import jctools
import pandas as pd
from Bio import SeqIO
from Bio import PDB
import os
import argparse

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
def read_noHead_data(inFile):
    dfRef_ext = inFile.split('/')[-1].split('.')[-1].lower()
    if dfRef_ext == 'csv':
        dfRef = pd.read_csv(inFile, comment='#', dtype={'PDB': object}, header=None)
    else:
        dfRef = pd.read_csv(inFile, sep='\t', comment='#', dtype={'PDB': object}, header=None)
    return dfRef

#------------------------------------------------------------------------------
# Parameters
# python PDBChains.py -i 14gs
# -o /home/jj140/scratch/ppi/batch/EC2/outputs
# -r /home/jj140/scratch/ppi/dbs/pdb_chain_enzyme.csv
# -s 10078
# -e 40899
#------------------------------------------------------------------------------
db_cath = '/Users/jjeong/local/project_dev/ppi/outputs/cath/cath_pdb.tsv'
db_ec = '/Users/jjeong/local/project_dev/ppi/outputs/cath/ec2_pdb.tsv'
db_go = '/Users/jjeong/local/project_dev/ppi/outputs/cath/go_pdb.tsv'
db_pfam = '/Users/jjeong/local/project_dev/ppi/outputs/cath/pfam_pdb.tsv'
db_scop = '/Users/jjeong/local/project_dev/ppi/outputs/cath/scop_pdb.tsv'


#------------------------------------------------------------------------------
# Search PDB and extract corresponding Chains
#------------------------------------------------------------------------------
dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})

#-- get file extension
dfRef_ext = pdbFile.split('/')[-1].split('.')[-1].tolower
if dfRef_ext == 'csv':
    dfRef = pd.read_csv(refFile, comment='#', dtype={'PDB': object})
else:
    dfRef = pd.read_csv(refFile, sep='\t', comment='#', dtype={'PDB': object})

#refFile='/home/jj140/scratch/ppi/dbs/pdb_chain_go.tsv'
#pdb_id='2c2k'
#Chain = pd.read_csv(refFile, sep='\t', comment='#', dtype={'PDB': object})
#for dfRef in Chunks:
#    Chains = ''.join(dfRef.loc[dfRef['PDB'] == pdb_id]['CHAIN'].unique().tolist())
#    print(Chains)

    
selPDB = dfPDB['PDB'][s_idx:e_idx].tolist()

for pdb_id in selPDB:
    Chains = ''.join(dfRef.loc[dfRef['PDB'] == pdb_id]['CHAIN'].unique().tolist())
    outFile = '{}/{}.tsv'.format(outDir, pdb_id)
    with open(outFile, 'w') as f:
        line = '{}\n'.format('\t'.join([pdb_id, Chains]))
        f.write(line)


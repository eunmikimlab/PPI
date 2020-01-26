"""
   ::Purpose::
   Extract PDB Chains corresponding to PDB name and write them into a file
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
# python PDBChains.py -i 14gs
# -o /home/jj140/scratch/ppi/batch/EC2/outputs
# -r /home/jj140/scratch/ppi/dbs/pdb_chain_enzyme.csv
# -s 10078
# -e 40899
#------------------------------------------------------------------------------
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help='A file list PDB ID that can be searched', action='store', required=True)
parser.add_argument("-o", "--outdir", help='Output directory to store PDB and Chain pairs', action='store', required=True)
parser.add_argument("-r", "--ref", help='Reference file where the PDB ID is searched', action='store', required=True)
parser.add_argument("-s", "--start", help='The index of PDB data to be started', action='store', required=True)
parser.add_argument("-e", "--end", help='The index of PDB data to be ended', action='store', required=True)

args = parser.parse_args()

pdbFile = args.input
outDir = args.outdir
refFile = args.ref
s_idx = int(args.start)
e_idx = int(args.end)

"""
outDir = '/Users/jjeong/local/project_dev/ppi/slurm/EC'
pdbFile = '/Users/jjeong/local/project_dev/ppi/dbs/outputs/unique_human_pdb.csv'
refFile = '/Users/jjeong/local/project_dev/ppi/dbs/pdb_chain_enzyme.csv'
s_idx = int(0)
e_idx = int(100)
"""

#------------------------------------------------------------------------------
# Search PDB and extract corresponding Chains
#------------------------------------------------------------------------------
dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})

#-- get file extension
dfRef_ext = refFile.split('/')[-1].split('.')[-1].lower()
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
    strChain = map(str,dfRef.loc[dfRef['PDB'] == pdb_id]['CHAIN'].unique().tolist())
    Chains = ''.join(strChain)
    outFile = '{}/{}.tsv'.format(outDir, pdb_id)
    with open(outFile, 'w') as f:
        line = '{}\n'.format('\t'.join([pdb_id, Chains]))
        f.write(line)


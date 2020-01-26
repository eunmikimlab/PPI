"""
   ::Purpose::
   Write lists of dssp and hssp download list to use rsync
   rsync -auvh --progress --files-from=/path/to/files.txt / /destination/path/
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 20, 2018"


import jctools
import pandas as pd
from Bio import SeqIO
from Bio import PDB
import os
import argparse

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
dbName = 'EC'
prefix = 'rsync://rsync.cmbi.ru.nl/dssp'
pdbFile = '/home/jj140/scratch/ppi/batch/{}_mergedDBs.csv'.format(dbName)   # PDB and all associated chains
out_dssp = '/home/jj140/scratch/ppi/batch/dsspList.txt'
out_hssp = '/home/jj140/scratch/ppi/batch/hsspList.txt'

dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})
selPDB = dfPDB['PDB'].tolist()

with open(out_dssp, 'w') as f:
    for pdb_id in selPDB:
        dssp = '{}/{}.dssp\n'.format(prefix, pdb_id.lower())
        f.write(dssp)
        
with open(out_hssp, 'w') as f:
    for pdb_id in selPDB:
        hssp = '{}/{}.hssp.bz2\n'.format(prefix, pdb_id.lower())
        f.write(hssp)

    
    
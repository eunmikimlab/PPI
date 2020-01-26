"""
   ::Purpose::
   This codes Extract PDBs that are corresponding to HUMAN
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 12, 2018"


import jctools
import pandas as pd
from Bio import SeqIO
from Bio import PDB
import os
#------------------------------------------------------------------------------
# Define Parameters
#------------------------------------------------------------------------------
DB_HOME = '/Users/jjeong/local/project_dev/ppi/dbs'
Texa = '{}/pdb_chain_taxonomy.tsv'.format(DB_HOME)
Data = '{}/pdb_chain_enzyme.csv'.format(DB_HOME)

OUT_HOME = '{}/outputs'.format(DB_HOME)
try:
    os.makedirs(OUT_HOME)
except OSError:
    pass

dfTexa = pd.read_csv(Texa, sep='\t', comment='#', dtype={'PDB': object, 'TAX_ID': int})
dfData = pd.read_csv(Data, comment='#', dtype={'PDB': object})

#-- Extract Human Only
dfHuman = dfTexa.loc[dfTexa['TAX_ID'] == 9606]
out_human = '{}/taxonomy_human_only.csv'.format(OUT_HOME)
dfHuman.to_csv(out_human, index=False)

#-- Write a file 'PDB_ID & Chains'
pdbIDs = dfHuman['PDB'].unique().tolist()
dfIDs = pd.DataFrame(pdbIDs, columns=['PDB'])
out_pdbid = '{}/unique_human_pdb.csv'.format(OUT_HOME)
dfIDs.to_csv(out_pdbid, index=False)


#uqPDBs = [[pdb, ''.join(dfHuman.loc[dfHuman['PDB'] == pdb]['CHAIN'].unique().tolist())] for pdb in pdbIDs]

#for pdb in pdbIDs:
#    tmp = dfHuman.loc[dfHuman['PDB'] == pdb]
#    tmpChain = tmp['CHAIN'].unique().tolist()
#    uqPDBs.append([pdb, ''.join(tmpChain)])

#dfPDBs = pd.DataFrame(uqPDBs, columns=['PDB', 'CHAIN'])
#dfPDBs.to_csv('PDB_CHAINS.csv', index=False)

print(Texa)
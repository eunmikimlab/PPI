"""
   ::Purpose::
   Inner Join PDB and Chain pairs among DBs
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 13, 2018"


import pandas as pd

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
cath ='/home/jj140/scratch/ppi/batch/cath/cath_pdb.tsv'
ec = '/home/jj140/scratch/ppi/batch/ec/ec_pdb.tsv'
go = '/home/jj140/scratch/ppi/batch/go/go_pdb.tsv'
pfam = '/home/jj140/scratch/ppi/batch/pfam/pfam_pdb.tsv'
scop = '/home/jj140/scratch/ppi/batch/scop/scop_pdb.tsv'
dbName = 'SCOP'
output = '/home/jj140/scratch/ppi/batch/{}_mergedDBs.csv'.format(dbName)
min_ch = 2      # minimum number of chains in PDB
max_ch = 3      # maximum number of chains in PDB

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
dfCath = pd.read_csv(cath, sep='\t', comment='#', header=None)
dfCath.columns = ['PDB', 'CHAIN']

dfEc = pd.read_csv(ec, sep='\t', comment='#', header=None)
dfEc.columns = ['PDB', 'CHAIN']

dfGo = pd.read_csv(go, sep='\t', comment='#', header=None)
dfGo.columns = ['PDB', 'CHAIN']

dfPfam = pd.read_csv(pfam, sep='\t', comment='#', header=None)
dfPfam.columns = ['PDB', 'CHAIN']

dfScop = pd.read_csv(scop, sep='\t', comment='#', header=None)
dfScop.columns = ['PDB', 'CHAIN']

dfDBs = [dfCath, dfEc, dfGo, dfPfam, dfScop]

#-- Do INNER JOIN based on PDB ID
dfMerged = reduce(lambda left, right: pd.merge(left, right, on=['PDB'], how='inner'), dfDBs)

dfMerged.columns = ['PDB', 'CATH', 'EC', 'GO', 'PFAM', 'SCOP']
noNaMerged = dfMerged.dropna()
Idx = (noNaMerged[dbName].str.len() >= min_ch) & (noNaMerged[dbName].str.len() <= max_ch)
selData = noNaMerged.loc[Idx]


selData.to_csv(output, index=False)

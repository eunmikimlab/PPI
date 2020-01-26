"""
   ::Purpose::
   To draw Heatmap, read data (e.g., out of 'searchMotif_seq.py') and add DBs family or group information
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 24, 2018"
__update__      = "Sep 25, 2018"


import pandas as pd
import os

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
dbName = 'scop'
selBase = 'seqBase'              # seqBase=sequence-based interface residues, strcBase=structure-based interface residues
selIntf = 'NoIntf'                 # target to search motif: 'Intf'= interface residue, 'NoIntf'=noninterface residues
len_motif = 3

#-- data output of 'searchMotif_seq.py'
if selBase == 'seqBase':
    inFile = '/Users/jjeong/local/project_dev/ppi/outputs/{}/seqBase/{}/{}_seqFreq_lenMotif_{}.csv'.format(dbName, selIntf, dbName.upper(), len_motif)
    OUT_HOME ='/Users/jjeong/local/project_dev/ppi/outputs/heatmap/seqBase/{}'.format(selIntf)

elif selBase == 'strcBase':
    inFile = '/Users/jjeong/local/project_dev/ppi/outputs/{}/strcBase/{}/{}_strcFreq_lenMotif_{}.csv'.format(dbName, selIntf, dbName.upper(), len_motif)
    OUT_HOME ='/Users/jjeong/local/project_dev/ppi/outputs/heatmap/strcBase/{}'.format(selIntf)


try:
    os.makedirs(OUT_HOME)
except OSError:
    pass

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------
dfInputs = pd.read_csv(inFile, comment='#', dtype={'PDB': object})
inIDs = dfInputs['CHAIN'].tolist()

if dbName.upper() == 'EC':
    dbParse = '/Users/jjeong/local/project_dev/ppi/dbs/pdb_chain_enzyme.csv'
    dfParse = pd.read_csv(dbParse, comment='#', dtype={'PDB': object})
    psIDs = dfParse['PDB']

    Category = []
    cnt = 0
    for id in inIDs:
        pdb = id.split('_')[0]
        chain = id.split('_')[2].split('.')[0]
        tmp = dfParse.loc[(dfParse['PDB'] == pdb) & (dfParse['CHAIN'] == chain[0])]['EC_NUMBER'].tolist()[0]
        Category.append(tmp)

        # progress
        cnt = cnt + 1
        if cnt % 100 == 0:
            print('[{}/{}] {} with len {} has been completed!'.format(cnt, len(inIDs), dbName, len_motif))

elif dbName.upper() == 'SCOP':
    dbParse = '/Users/jjeong/local/project_dev/ppi/dbs/pdb_chain_scop_uniprot.csv'
    dfParse = pd.read_csv(dbParse, comment='#', dtype={'PDB': object})
    psIDs = dfParse['PDB']

    dbParse2 = '/Users/jjeong/local/project_dev/ppi/dbs/parse/dir.des.scop_1.75.txt'
    dfParse2 = pd.read_csv(dbParse2, comment='#', sep='\t', header=None, names=['SUNID', 'DTYPE', 'SCOP_ID', 'DOMAIN', 'INFO'])

    Category = []
    cnt = 0
    for id in inIDs:
        pdb = id.split('_')[0]
        chain = id.split('_')[2].split('.')[0]
        sunid = dfParse.loc[(dfParse['PDB'] == pdb) & (dfParse['CHAIN'] == chain[0])]['SUNID'].tolist()[0]
        tmp = dfParse2.loc[(dfParse2['SUNID'] == sunid)]['SCOP_ID'].tolist()[0]
        Category.append(tmp)

        # progress
        cnt = cnt + 1
        if cnt % 100 == 0:
            print('[{}/{}] {} with len {} has been completed!'.format(cnt, len(inIDs), dbName, len_motif))


#-- Store Category Information
dfInputs['CATEGORY'] = Category

if selBase == 'seqBase':
    outFile = '{}/heatmap_{}_seqFreq_lenMotif_{}.csv'.format(OUT_HOME, dbName.upper(), len_motif)

elif selBase == 'strcBase':
    outFile = '{}/heatmap_{}_strcFreq_lenMotif_{}.csv'.format(OUT_HOME, dbName.upper(), len_motif)

dfInputs.to_csv(outFile, index=False)

print('== DONE ==')

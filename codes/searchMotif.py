"""
   ::Purpose::
   Search sequence and returns matched index
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 21, 2018"


import jctools
import pandas as pd
import collections
import itertools
import os
import re
import math

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)


#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
dbName = 'ec'
len_motiff = 2
lbl_fext = '.csv$'
LABEL_HOME = '/Users/jjeong/local/project_dev/ppi/outputs/labels'
pdbFile = '/Users/jjeong/local/project_dev/ppi/outputs/{}_mergedDBs.csv'.format(dbName.upper())              # Unique PDB list obtained from 'joinDBs.py'

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

#-- 20 AAs
AAs=['V', 'L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', 'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E', 'N', 'D']

#-- calculating Motif. Motif is all possible combinations of AAs with length 'len_motiff'
if len_motiff < 2:
    Motif = sorted(AAs)
else:
    Motif = list(itertools.combinations_with_replacement(AAs, len_motiff))
    Motif = sorted(list(set(map(lambda x: ''.join(sorted(list(x))), Motif))))

#-- Search all label files
tg_dir = '{}/{}/outputs'.format(LABEL_HOME, dbName)
lblFiles = [f for f in os.listdir(tg_dir) if (os.path.isfile(os.path.join(tg_dir, f)) & (re.search(lbl_fext, f) is not None))]

#-- Read Chain Information
dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})

#-- Read and count
dbName = dbName.upper()
outData = []
cnt = 0
for l_file in lblFiles:
    cnt = cnt + 1
    pdb_id = l_file.split('_')[0]
    Chains = ''.join(dfPDB.loc[dfPDB['PDB'] == pdb_id][dbName].unique().tolist())
    num_chain = len(Chains)
    pdb = '{}/{}'.format(tg_dir, l_file)
    Label = pd.read_csv(pdb)

    #-- extract interface residues only
    #Seq = ''.join(Label['SEQ'].tolist())
    Seq = ''.join(Label.loc[Label['CLASS'] == 1]['SEQ'].tolist())

    # count motiff
    motif, freq = jctools.countMotif(Seq, Motif)

    # normalize with sequence Length and nC2 n=number of chains, 2= always considering pairs
    if len(Seq) > 0:
        normF = map(lambda x: x / float(len(Seq)) * float(nCr(len(Chains), 2)), freq)
    else:
        normF = [0.0] * len(Motif)
    outData.append(normF)

    if cnt % 100 == 0:
        print('[{}/{}] has been completed!'.format(cnt, len(lblFiles)))
        break

dfData = pd.DataFrame(outData, columns=motif, dtype=float)

print(pdb)
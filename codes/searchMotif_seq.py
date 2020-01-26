"""
   ::Purpose::
   Calculated normalized Motif frequence within interface residues in Sequence Order
   Sequence Order calculates Motif frequence based on individual interface residues
   E.g., NNNNIIINNN : N=non interface residues, I=interfarce residues
   With Motif length 3
   'NNI', 'NII', 'INN' will be considered to the motif frequency of 5th 'I'
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 21, 2018"
__update__      = "Sep 23, 2018"



import jctools
import pandas as pd
import collections
import itertools
import os
import re
import math
import numpy as np
#import matplotlib.pyplot as plt

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
dbName = 'go'     # ec, scop, pfam, cath, go
len_motif = 4
sel_class = 0       # sel_class: 0=noninterface residues, 1=interface residues
lbl_fext = '.csv$'
LABEL_HOME = '/Users/jjeong/local/project_dev/ppi/outputs/labels'
pdbFile = '/Users/jjeong/local/project_dev/ppi/outputs/{}_mergedDBs.csv'.format(dbName.upper())              # Unique PDB list obtained from 'joinDBs.py'

if sel_class == 0:
    OUT_HOME = '/Users/jjeong/local/project_dev/ppi/outputs/{}/seqBase/NoIntf'.format(dbName)
elif sel_class == 1:
    OUT_HOME = '/Users/jjeong/local/project_dev/ppi/outputs/{}/seqBase/Intf'.format(dbName)

try:
    os.makedirs(OUT_HOME)
except OSError:
    pass

#------------------------------------------------------------------------------
# Functions
#------------------------------------------------------------------------------
def nCr(n,r):
    f = math.factorial
    return f(n) / f(r) / f(n-r)

# Count frequency of Motif by searching before and after the position of interface residues
# Label must include 'CLASS' and 'SEQ' columns with row index start with 0
# e.g., Label =
#    CHAIN  INDEX SEQ  CLASS
# 0       A      0   S      0
# 1       A      1   M      0
# 2       A      2   E      1
# 3       A      3   N      0
#----------------------------------------------------------
def countSeqIRs(Label, Motif):
    len_motif = len(Motif[0])
    s_Ridx = Label.index.tolist()[0]  # first residue index
    e_Ridx = Label.index.tolist()[-1] # last residue index
    dfIRs = Label.loc[Label['CLASS'] == 1]['SEQ']
    Idx = dfIRs.index.tolist()
    listFreq = []
    for i in Idx:
        if (i - len_motif <s_Ridx) and (i + len_motif - 1 <= e_Ridx):
            st_i = s_Ridx
            ed_i = i + len_motif - 1

        elif (i - len_motif < s_Ridx) and (i + len_motif - 1 > e_Ridx):
            st_i = s_Ridx
            ed_i = e_Ridx

        elif (i - len_motif >= s_Ridx) and (i + len_motif - 1 > e_Ridx):
            st_i = i - len_motif + 1
            ed_i = e_Ridx
        else:
            st_i = i - len_motif + 1
            ed_i = i + len_motif - 1

        tgSeq = ''.join(Label.iloc[st_i:ed_i]['SEQ'].tolist())

        # count motiff
        motif, freq = jctools.countMotif(tgSeq, Motif)

        # Append frequency
        listFreq.append(freq)

    # sum all frequencies
    totalFreq = map(sum, zip(*listFreq))
    return totalFreq


#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

#-- 20 AAs
AAs=['V', 'L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', 'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E', 'N', 'D']

#-- calculating Motif. Motif is all possible combinations of AAs with length 'len_motiff'
if len_motif < 2:
    Motif = sorted(AAs)
else:
    Motif = list(itertools.combinations_with_replacement(AAs, len_motif))
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
PDBID = []
CHAIN = []
for l_file in lblFiles:
    cnt = cnt + 1
    pdb_id = l_file.split('_')[0]
    Chains = ''.join(dfPDB.loc[dfPDB['PDB'] == pdb_id][dbName].unique().tolist())

    PDBID.append(pdb_id)
    CHAIN.append(l_file)

    num_chain = len(Chains)
    pdb = '{}/{}'.format(tg_dir, l_file)
    Label = pd.read_csv(pdb)

    # extract interface residues only CLASS = 1
    #extract NON interface residues only CLASS = 0
    #Seq = ''.join(Label['SEQ'].tolist())
    Seq = ''.join(Label.loc[Label['CLASS'] == sel_class]['SEQ'].tolist())

    # count motif
    # motif, freq = jctools.countMotif(Seq, Motif)
    if len_motif > 1:
        freq = countSeqIRs(Label, Motif)
    else:
        motif, freq = jctools.countMotif(Seq, Motif)

    # normalize with sequence Length and nC2 n=number of chains, 2= always considering pairs
    if len(Seq) > 0:
        normF = map(lambda x: x / float(len(Seq)) * float(nCr(len(Chains), 2)), freq)
    else:
        normF = [0.0] * len(Motif)
    outData.append(normF)

    if cnt % 100 == 0:
        print('[{}/{}] {} with len {} has been completed!'.format(cnt, len(lblFiles), dbName, len_motif))


dfData = pd.DataFrame(outData, columns=Motif, dtype=float)
dfData['PDB'] = PDBID
dfData['CHAIN'] = CHAIN
outFile = '{}/{}_seqFreq_lenMotif_{}.csv'.format(OUT_HOME, dbName, len_motif)
dfData.to_csv(outFile, index=False)

#plt.pcolor(dfData)
#plt.yticks(np.arange(0.5, len(dfData.index), 1), dfData.index)
#plt.xticks(np.arange(0.5, len(dfData.columns), 1), dfData.columns)
#plt.show()

print ('== DONE ==')

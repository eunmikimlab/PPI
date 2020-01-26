import jctools
import pprint
import pandas as pd
from Bio import SeqIO
from Bio import PDB

import itertools

import re
from StringIO import StringIO
import argparse
import json
import requests
import time


pp = pprint.PrettyPrinter(indent=4)

#------------------------------------------------------------------------------
# Basic usage of windows classes in 'jc_window.py'
#------------------------------------------------------------------------------
mySeq = '1234567890'
objWin = jctools.winSeq(mySeq)

myWin = objWin.linearwin(5, 2, 0)
SEQs = jctools.retrieveSeq(myWin, len(mySeq))

myWin2 = objWin.periodicWin(21, '-1')
SEQ2 = jctools.retrieveSeq(myWin2, len(mySeq))


#-- convert sequence to mapped property
finput = '/Users/jjeong/local/project_dev/ppi/inputs/AAprofiles.csv'
ifasta ='/Users/jjeong/local/project_dev/ppi/inputs/uniprot_human.fasta'
ipdb   = '/Users/jjeong/local/project_dev/ppi/inputs/pdb/1fvc.pdb'
aaRef  = '/Users/jjeong/local/project_dev/ppi/inputs/AAprofiles.csv'

#-- read AAprofiles
df = pd.read_csv(finput)

#-- Read FASTA
with open(ifasta, 'rU') as handle:
    for record in SeqIO.parse(handle, 'fasta'):
        print("Record ID = %s " % (record.id))
        print("Sequence  = %s " % (record.seq))
        print("Alphabet  = %s " % (record.seq.alphabet))
        break


#-- Read PDB : R3 = 3 letter AA name, R1= 1 letter AA name
R3, R1 = jctools.pdb_read_seq(ipdb)

"""
with open(ipdb, 'rU') as handle:
    for record in SeqIO.parse(handle, 'pdb-seqres'):
        print("Record ID = %s " % (record.id))
        print("Chain ID = %s" % (record.annotations['chain']))
        print("Sequence  = %s " % (record.seq))
        print("Alphabet  = %s " % (record.seq.alphabet))
"""



#------------------------------------------------------------------------------
# Read sequence from a file and slide window & map the physicochemical values
#------------------------------------------------------------------------------

#-- Mapp sequence with physicochemical properties
#aaSeq = record.seq                          # get AA sequence
aaSeq = R1['A']


# sliding window on a sequence
objWin = jctools.winSeq(aaSeq)
pdWin = objWin.periodicWin(5, '-1')                 # fill outside heading and tailing window with '-1'
pdSEQ = jctools.retrieveSeq(pdWin, len(aaSeq))      # get sequence as a list

#-- read reference sequence (physicochemcial properties)
ref = pd.read_csv(aaRef)

#-- convert slide sequence with designated properity (e.g., 'MS', 'HPO')
# Property name is case sensetive, so MUST match the COLUMN HEAD in './input/AAprofiles.csv'
newSEQ = []
for w in pdSEQ:
    newSEQ.append(jctools.convertAA(w, ref, 'MS'))


mapData = pd.DataFrame(newSEQ)


mySeq = 'EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS'





#------------------------------------------------------------------------------
# Run XSSP to get HSSP results from a protein sequence
#------------------------------------------------------------------------------
"""
#- using AA sequence
hssp = '/Users/jjeong/local/project_dev/ppi/outputs/hssp_test2.txt'
hsspStr = jctools.seq_to_hssp(mySeq)
with open(hssp, 'w') as f:
    f.write(hsspStr)
"""
#- using PDB ID
hssp = '/Users/jjeong/local/project_dev/ppi/outputs/hssp_5v9o.txt'
hsspStr = jctools.pdbid_to_hssp('1L7I')
with open(hssp, 'w') as f:
    f.write(hsspStr)

#------------------------------------------------------------------------------
# Read HSSP file and write Feature Matrix
#------------------------------------------------------------------------------
hssp = '/Users/jjeong/local/project_dev/ppi/outputs/hssp_mySeq.txt'
hsspProfile = jctools.profile_hssp(hssp)

# get all column headers
#hsspProfile.columns.tolist()

# get AA order
aaOrder = hsspProfile.columns.tolist()[2:22]

Data = jctools.build_fMtx_hssp(hssp=hsspProfile, seq_len=len(mySeq), win_size=21)

# write data
out_path = '/Users/jjeong/local/project_dev/ppi/outputs/hssp_feature_matrix.csv'
Data.to_csv(out_path)


#------------------------------------------------------------------------------
# Map Seq to Physicochemical properties and build/write Feature Matrix with average
#------------------------------------------------------------------------------
aaRef  = '/Users/jjeong/local/project_dev/ppi/inputs/AAprofiles.csv'

tmpAve = jctools.build_fMtx_refAve(Seq=mySeq, AAref=aaRef, Attribute='HPO', win_size=11)
tmp = jctools.build_fMtx_Attr(Seq=mySeq, AAref=aaRef, Attribute='HPO', win_size=11)
print(tmp)


#------------------------------------------------------------------------------
# Map Seq to Physicochemical properties and build/write Feature Matrix with momentum
#------------------------------------------------------------------------------
mtAve = jctools.build_fMtx_refmtAve(Seq=mySeq, AAref=aaRef, Attribute='MS')

#------------------------------------------------------------------------------
# Create assign vector to keep the sequence information
#------------------------------------------------------------------------------
seqMtx = jctools.build_fMtx_Seq(mySeq)


#------------------------------------------------------------------------------
# Search the 20 AAs distance from left and right of the current residue
# None = -1, Myself = 0
#------------------------------------------------------------------------------

#-- using whole sequence
winLeft, winRight = jctools.build_fMtx_SearchDist(mySeq)

#-- using windowing
winLeft2, winRight2 = jctools.build_fMtx_SearchDistWin(mySeq, win_size=11)


#------------------------------------------------------------------------------
# Retrieve PDB and store it into a local driver
# https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
#------------------------------------------------------------------------------
pdb_id = '6DM0'
PDB_HOME = '/Users/jjeong/local/project_dev/ppi/outputs/pdb'
pdbDL = PDB.PDBList()
pdbDL.retrieve_pdb_file(pdb_code=pdb_id, pdir=PDB_HOME)

# read any file format
cif_file = "{}/{}.cif".format(PDB_HOME, pdb_id.lower())
cifParser = PDB.MMCIFParser()
structure = cifParser.get_structure(pdb_id, cif_file)

pdb_file = "{}/{}.pdb".format(PDB_HOME, pdb_id.lower())
pdbIO = PDB.PDBIO()
pdbIO.set_structure(structure[0])
pdbIO.save(pdb_file)

#------------------------------------------------------------------------------
# Get Surface residue
# https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ
# https://biopython.org/DIST/docs/api/Bio.PDB.DSSP%27-module.html
# Download and Install DSSP is required
# ftp://ftp.cmbi.ru.nl//pub/molbio/software/dssp-2/
#------------------------------------------------------------------------------
# Read PDB
#pdb_file = '/Users/jjeong/local/project_dev/ppi/codes/utils/1A2Y.pdb'

# load structure
pdbParser = PDB.PDBParser()
structure = pdbParser.get_structure(pdb_id, pdb_file)
model = structure[0]
#dssp = PDB.DSSP(model=model, in_file="/Users/jjeong/local/project_dev/ppi/outputs/pdb/6dm0.dssp", file_type='DSSP')
dssp = PDB.DSSP(model=model, in_file=pdb_file, dssp='mkdssp', acc_array='Sander', file_type='PDB')

#-- To see Max ACC
maxAcc = dssp.residue_max_acc



print(hsspProfile)


"""
   ::Purpose::
   This package is implemented for generating Amino Acids (AA) features by read input file (i.e., ./inputs/AAprofiles.csv)
   Required packages
   1. pandas
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Feb 13, 2018"


import pandas as pd
import re
from StringIO import StringIO
import json
import requests
import time
import jc_window
import numpy as np
from Bio import PDB
import ntpath


#----------------------------------------------------------
# Get file name from the path
# This will get the last item in the path evenif it has '/' in file name
#----------------------------------------------------------
def get_filename(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

#----------------------------------------------------------
# Read protein sequence from PDB
#----------------------------------------------------------
def pdb_read_seq(pdb_file=None, model_idx=0):
    Ref = pd.DataFrame({"LETTER1": ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'],
                        "LETTER3": ['ALA', 'CYS', 'ASP', 'GLU', 'PHE', 'GLY', 'HIS', 'ILE', 'LYS', 'LEU', 'MET', 'ASN', 'PRO', 'GLN', 'ARG', 'SER', 'THR', 'VAL', 'TRP', 'TYR']})

    pdb_id = get_filename(pdb_file).split('.')[0]

    pdbParser = PDB.PDBParser()
    structure = pdbParser.get_structure(pdb_id, pdb_file)
    model = structure[model_idx]

    #-- Get list of Chains
    Chains = [c._id for c in model.get_chains()]

    #-- Read Sequence from model
    Res3 = {}
    for c in Chains:
        tmp = []
        for residue in model[c].get_residues():
            tmp.append((residue.get_resname(), residue.id[1]))
        Res3[c] = tmp


    Res1 = {}
    for c in Chains:
        tmp = []
        cnt = 0
        for r3 in Res3[c]:
            cnt = cnt + 1
            #-- pdb can contain more than 20 AAs (e.g., CYZ - ligand)
            if r3[0] in Ref['LETTER3'].values:
                tmp.append((str(Ref[Ref['LETTER3'] == r3[0]]['LETTER1'].item()), r3[1]))
            else:
                print("<!!Warning!!> In chain '{}', '{}:{}' is not the one of 20 amino acids.".format(c, r3[0], cnt))

        Res1[c] = tmp

    return Res3, Res1

#----------------------------------------------------------
# Mapping Amino Acid Residues to Physicochemical properties
# 1) Seq is an array of AA sequence: list('ACDEG')
# 2) Ref is an PANDAS object of CSV file (i.e., ./input/AAprofiles.csv)
# 3) Pro is the name of Property that is converted to a sequence
# 4) dummy is a dummy character that will be ignored
#   (e.g., caused by leading and tailing of windows on a sequence)
# 5) aaName is
#----------------------------------------------------------
def convertAA(Seq, Ref, Prop='HPO', dummy='-1', aaName='LETTER1'):
    newSEQ = []
    for aa in Seq:
        if aa == dummy:
            tmp = str(dummy)
        else:
            tmp = str(Ref[Ref[aaName] == aa.upper()][Prop].item())
        newSEQ.append(tmp)
    return(newSEQ)


#------------------------------------------------------------------------------
# Using HSSP
#url_create = 'http://www.cmbi.umcn.nl/xssp/api/create/sequence/hssp_hssp/'
#------------------------------------------------------------------------------
"""
# HSSP from Amino Acid Sequences
"""
def seq_to_hssp(Seq, rest_url='http://www.cmbi.umcn.nl/xssp'):
    #Seq=mySeq
    #rest_url='http://www.cmbi.umcn.nl/xssp'

    api_create ='api/create/sequence/hssp_hssp'
    api_status ='api/status/sequence/hssp_hssp'
    api_result ='api/result/sequence/hssp_hssp'

    url_create = '{}/{}/'.format(rest_url, api_create)
    url_status = '{}/{}/'.format(rest_url, api_status)
    url_result = '{}/{}/'.format(rest_url, api_result)

    #-- submit job
    r = requests.post(url_create, data={'data':Seq})
    r.raise_for_status()
    job_id = json.loads(r.text)['id']
    print "Job submitted successfully. Id is: '{}'".format(job_id)

    #-- loop until results are ready
    ready = False
    while not ready:
        ask_status = '{}/{}/'.format(url_status, job_id)
        r = requests.get(ask_status)
        r.raise_for_status()

        #-- check status
        status = json.loads(r.text)['status']
        print "Job status is: '{}'".format(status)
        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE', 'REVOKED']:
            raise Exception(json.loads(r.text)['message'])
        else:
            time.sleep(5)
    else:
        ask_result = '{}/{}/'.format(url_result, job_id)
        r = requests.get(ask_result)
        r.raise_for_status()
        result = json.loads(r.text)['result']

        # Return the result to the caller, which prints it to the screen.
        return result



"""
# HSSP from PDB file inputs
"""
def pdbfile_to_hssp(pdb_file_path, rest_url='http://www.cmbi.umcn.nl/xssp'):
    # Read the pdb file data into a variable
    files = {'file_': open(pdb_file_path, 'rb')}

    api_create ='api/create/pdb_file/hssp_hssp'
    api_status ='api/status/pdb_file/hssp_hssp'
    api_result ='api/result/pdb_file/hssp_hssp'

    url_create = '{}/{}/'.format(rest_url, api_create)
    url_status = '{}/{}/'.format(rest_url, api_status)
    url_result = '{}/{}/'.format(rest_url, api_result)

    #url_create = '{}api/create/pdb_file/hssp_hssp/'.format(rest_url)
    r = requests.post(url_create, data={'data':Seq})
    r.raise_for_status()
    job_id = json.loads(r.text)['id']
    print "Job submitted successfully. Id is: '{}'".format(job_id)

    # Loop until the job running on the server has finished, either successfully
    # or due to an error.
    ready = False
    while not ready:
        ask_status = '{}/{}/'.format(url_status, job_id)
        r = requests.get(ask_status)
        r.raise_for_status()

        status = json.loads(r.text)['status']
        print "Job status is: '{}'".format(status)

        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE', 'REVOKED']:
            raise Exception(json.loads(r.text)['message'])
        else:
            time.sleep(5)
    else:
        ask_result = '{}/{}/'.format(url_result, job_id)
        r = requests.get(ask_result)
        r.raise_for_status()
        result = json.loads(r.text)['result']

        # Return the result to the caller, which prints it to the screen.
        return result



"""
# HSSP from PDB ID
"""
def pdbid_to_hssp(pdb_id, rest_url='http://www.cmbi.umcn.nl/xssp'):
    # Read the pdb file data into a variable

    api_create ='api/create/pdb_id/hssp_hssp'
    api_status ='api/status/pdb_id/hssp_hssp'
    api_result ='api/result/pdb_id/hssp_hssp'

    url_create = '{}/{}/'.format(rest_url, api_create)
    url_status = '{}/{}/'.format(rest_url, api_status)
    url_result = '{}/{}/'.format(rest_url, api_result)

    r = requests.post(url_create, data={'data':pdb_id})
    r.raise_for_status()
    job_id = json.loads(r.text)['id']
    print "Job submitted successfully. Id is: '{}'".format(job_id)

    # Loop until the job running on the server has finished, either successfully
    # or due to an error.
    ready = False
    while not ready:
        ask_status = '{}/{}/'.format(url_status, job_id)
        r = requests.get(ask_status)
        r.raise_for_status()

        status = json.loads(r.text)['status']
        print "Job status is: '{}'".format(status)

        if status == 'SUCCESS':
            ready = True
        elif status in ['FAILURE', 'REVOKED']:
            raise Exception(json.loads(r.text)['message'])
        else:
            time.sleep(5)
    else:
        ask_result = '{}/{}/'.format(url_result, job_id)
        r = requests.get(ask_result)
        r.raise_for_status()
        result = json.loads(r.text)['result']

        # Return the result to the caller, which prints it to the screen.
        return result


#------------------------------------------------------------------------------
# Read HSSP file and extract hssp sequence profiles
# hssp   = path of HSSP file
# ptnStr = section that contains SEQUENCE PROFILE
#------------------------------------------------------------------------------
def profile_hssp(hssp, ptnStr = '## SEQUENCE PROFILE AND ENTROPY', colOrder=['SeqNo', 'PDBNo', 'V', 'L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', 'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E', 'N', 'D', 'NOCC', 'NDEL', 'NINS', 'ENTROPY', 'RELENT', 'WEIGHT']):
    buffStr = []
    keepFlag = False    # True if ptnStr is found
    with open(hssp, 'r') as f:
        for line in f:
            if keepFlag:
                #-- another comment '##' appears quit the loop
                if line.startswith('##'):
                    break
                else:
                    tmp = re.sub(' +', '\t', line.strip())
                    buffStr.append(tmp)

            # make keepFlag True and read next line and store the line into buffStr
            if line.startswith(ptnStr):
                keepFlag = True

    df = pd.read_table(StringIO('\n'.join(buffStr)))
    df = df[colOrder]
    return df


#------------------------------------------------------------------------------
# Build HSSP Data Matrix
# hssp = n x m where 'n' is AA sequence position and 'm' is 20 AAs + other hssp properties
#------------------------------------------------------------------------------
def build_fMtx_hssp(hssp, seq_len=None, win_size=21, null_val='-1', aa_start=2, aa_end=22):

    if seq_len == None :
        seq_len = hssp.shape[0]
    else:
        assert (hssp.shape[0] == seq_len), "The length of Sequence ({}) and the number of rows in HSSP ({}) are not matched".format(seq_len, hssp.shape[0])

    # get AA order
    #aaOrder = list(hsspProfile)[2:22]

    mySeqIdx = range(0, seq_len)
    objWin   = jc_window.winSeq(mySeqIdx)
    myWin    = objWin.periodicWin(win_size, null_val)
    seqIdx = jc_window.retrieveSeq(myWin, len(mySeqIdx))

    dataMtx = []
    for i in range(0, len(seqIdx)):
        winIdx = seqIdx[i]
        tmp = []
        for j in range(0, len(winIdx)):
            aaIdx = winIdx[j]
            if aaIdx == null_val:
                tmp = tmp + [float(null_val)]*(aa_end - aa_start)
            else:
                tmp = tmp + hssp.iloc[aaIdx, aa_start:aa_end].tolist()
        dataMtx.append(tmp)

    dfData = pd.DataFrame(dataMtx)

    return dfData



#------------------------------------------------------------------------------
# Attribute Average score per a window
#------------------------------------------------------------------------------
def build_fMtx_Attr(Seq=None, AAref=None, Attribute=None, win_size=11):
    assert (Seq != None), "Amino acid sequence must be given!"
    assert (AAref != None), "Reference table must be given!"
    assert (Attribute != None), "Attributes in Reference table must be given!"

    # sliding window on a sequence
    objWin = jc_window.winSeq(Seq)
    pdWin = objWin.periodicWin(win_size=win_size, pad='-1')  # fill outside heading and tailing window with '-1'
    pdSEQ = jc_window.retrieveSeq(pdWin, len(Seq))           # get sequence as a list

    # -- read reference sequence (physicochemcial properties)
    ref = pd.read_csv(AAref)

    # -- convert slide sequence with designated properity (e.g., 'MS', 'HPO')
    # Property name is case sensetive, so MUST match the COLUMN HEAD in './input/AAprofiles.csv'
    mapSEQ = []
    for w in pdSEQ:
        tmp = map(float, convertAA(w, ref, Attribute))
        mapSEQ.append(tmp)


    dfData = pd.DataFrame(mapSEQ)

    return dfData


#------------------------------------------------------------------------------
# Attribute Average score per a window
#------------------------------------------------------------------------------
def build_fMtx_refAve(Seq=None, AAref=None, Attribute=None, win_size=11, dummy='-999'):
    assert (Seq != None), "Amino acid sequence must be given!"
    assert (AAref != None), "Reference table must be given!"
    assert (Attribute != None), "Attributes in Reference table must be given!"

    # sliding window on a sequence
    objWin = jc_window.winSeq(Seq)
    pdWin = objWin.periodicWin(win_size=win_size, pad='-999')  # fill outside heading and tailing window with '-1'
    pdSEQ = jc_window.retrieveSeq(pdWin, len(Seq))           # get sequence as a list

    # -- read reference sequence (physicochemcial properties)
    ref = pd.read_csv(AAref)

    # -- convert slide sequence with designated properity (e.g., 'MS', 'HPO')
    # Property name is case sensetive, so MUST match the COLUMN HEAD in './input/AAprofiles.csv'
    mapSEQ = []
    for w in pdSEQ:
        tmp = map(float, convertAA(Seq=w, Ref=ref, Prop=Attribute, dummy=dummy))
        sel = [x for x in tmp if x != -999.0]
        ave = float(sum(sel)) / max(len(sel), 1)
        mapSEQ.append(ave)

    return mapSEQ



#------------------------------------------------------------------------------
# Attribute momentum Average score per a window
#------------------------------------------------------------------------------
def build_fMtx_refmtAve(Seq=None, AAref=None, Attribute=None, win_size=11, Angles=[60, 120, 180, 240, 300], dummy='-999'):
    assert (Seq != None), "Amino acid sequence must be given!"
    assert (AAref != None), "Reference table must be given!"
    assert (Attribute != None), "Attributes in Reference table must be given!"

    # sliding window on a sequence
    objWin = jc_window.winSeq(Seq)
    pdWin = objWin.periodicWin(win_size=win_size, pad='-999')  # fill outside heading and tailing window with '-1'
    pdSEQ = jc_window.retrieveSeq(pdWin, len(Seq))           # get sequence as a list

    #-- read reference sequence (physicochemcial properties)
    ref = pd.read_csv(AAref)

    #-- convert slide sequence with designated properity (e.g., 'MS', 'HPO')
    # Property name is case sensetive, so MUST match the COLUMN HEAD in './input/AAprofiles.csv'
    mapSEQ = []
    for w in pdSEQ:
        tmp = map(float, convertAA(Seq=w, Ref=ref, Prop=Attribute, dummy=dummy))
        sel = np.array([x for x in tmp if x != -999.0 ])
        m = np.zeros(len(sel))
        for i in Angles:
            for j in Angles:
                #m = m + (sel * np.power(np.sin(np.deg2rad(i)), 2) + sel * np.power(np.cos(np.deg2rad(j)), 2))
                m = m + np.power(sel * np.sin(np.deg2rad(i)), 2) + np.power(sel * np.cos(np.deg2rad(j)), 2)

        mave = m / max((i * j), 1)
        mave = sum(mave) / max(len(mave), 1)
        mapSEQ.append(mave.tolist())

    return mapSEQ



#------------------------------------------------------------------------------
# Attribute Average score per a window
#------------------------------------------------------------------------------
def build_fMtx_Seq(Seq=None, colOrder=['V', 'L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', 'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E', 'N', 'D']):
    Seq = Seq.upper()
    dataMtx = []
    for i in range(len(Seq)):
        tmp = [0] * len(colOrder)
        tmp[colOrder.index(Seq[i])] = 1
        dataMtx.append(tmp)
    dfData = pd.DataFrame(dataMtx)
    dfData.columns = colOrder

    return dfData


#------------------------------------------------------------------------------
#  Search the 20 AAs distance from LEFT and RIGHT of the current residue
#------------------------------------------------------------------------------
def build_fMtx_SearchDist(Seq=None, colOrder=['V', 'L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', 'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E', 'N', 'D']):
    mtxLeft = []
    mtxRight = []
    for i in range(len(Seq)):
        #-- search left side of current residue including itself
        lstLeft = []
        leftR = Seq[0:i+1]
        revR = leftR[::-1]
        for res in colOrder:
            try:
                idx = revR.index(res)
            except ValueError:
                idx = -1
            lstLeft.append(idx)
        mtxLeft.append(lstLeft)

        # -- search right side of current residue including itself
        lstRight = []
        rightR = Seq[i:len(Seq)]
        for res in colOrder:
            try:
                idx = rightR.index(res)
            except ValueError:
                idx = -1
            lstRight.append(idx)
        mtxRight.append(lstRight)

    dfLeft = pd.DataFrame(mtxLeft)
    dfRight = pd.DataFrame(mtxRight)

    dfLeft.columns = colOrder
    dfRight.columns = colOrder

    return dfLeft, dfRight


#------------------------------------------------------------------------------
#  Merging Pandas list of lists (2) by row
#------------------------------------------------------------------------------
def pd2MergeRow(pdData=None):
    tmp = []
    for i in range(len(pdData)):
        tmp.append(pdData.iloc[i, ].values.tolist())
    flatList = [item for sub1 in tmp for item in sub1 ]

    return flatList

#------------------------------------------------------------------------------
#  Search the 20 AAs distance from LEFT and RIGHT of the current residue
#------------------------------------------------------------------------------
def build_fMtx_SearchDistWin(Seq=None, win_size=11, colOrder=['V', 'L', 'I', 'M', 'F', 'W', 'Y', 'G', 'A', 'P', 'S', 'T', 'C', 'H', 'R', 'K', 'Q', 'E', 'N', 'D']):
    # -- window based distance
    objWin = jc_window.winSeq(Seq)
    pdWin = objWin.periodicWin(win_size=win_size, pad='-1')  # fill outside heading and tailing window with '-1'
    pdSEQ = jc_window.retrieveSeq(pdWin, len(Seq))  # get sequence as a list
    flatLeft = []
    flatRight = []
    for w in pdSEQ:
        winLeft, winRight = build_fMtx_SearchDist(w, colOrder=colOrder)
        flatLeft.append(pd2MergeRow(winLeft))
        flatRight.append(pd2MergeRow(winRight))

    winLeft = pd.DataFrame(flatLeft)
    winRight = pd.DataFrame(flatRight)

    return winLeft, winRight

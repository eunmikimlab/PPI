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
import os
import itertools
import copy
import operator
import collections

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
            #else:
            #    print("<!!Warning!!> In chain '{}', '{}:{}' is not the one of 20 amino acids.".format(c, r3[0], cnt))

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



###############################################################################
# Build Contact, Interface, and Surface Residue Labels
# Some of codes are from
# https://www.blopig.com/blog/2013/10/get-pdb-intermolecular-protein-contacts-and-interface-residues/
# http://www.stats.ox.ac.uk/~krawczyk/GetContacts.zip
###############################################################################

#--------------------------------------------------------------------
# Check numbers
#--------------------------------------------------------------------
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


#--------------------------------------------------------------------
# Get atom list with information of Chain, Resid, and icode
#--------------------------------------------------------------------
def get_atom_list(structure, chains):
    output = dict()
    for chain in structure:
        if chain.id in chains:
            for residue in chain.get_residues():
                hetflag, resseq, icode = residue.get_id()
                the_id = (chain.id + "_" + str(resseq) + "_" + icode).strip()
                for atom in residue.get_unpacked_list():
                    if hetflag == ' ':
                        if the_id in output:
                            output[the_id].append(atom)
                        else:
                            output[the_id] = [atom]
    return output


def get_filename(path):
    head, tail = ntpath.split(path)
    return tail or ntpath.basename(head)

#----------------------------------------------------------
# Read Chains from PDB
#----------------------------------------------------------
def pdb_read_chain(pdb_file=None, model_idx=0):
    pdb_id = get_filename(pdb_file).split('.')[0]

    pdbParser = PDB.PDBParser()
    structure = pdbParser.get_structure(pdb_id, pdb_file)
    model = structure[model_idx]

    #-- Get list of Chains
    Chains = [c._id for c in model.get_chains()]

    return Chains
    

def is_contact(res_1, other_atoms, cutoff):
    for atom in res_1:
        ns = PDB.NeighborSearch(other_atoms)
        center = atom.get_coord()
        neighbors = ns.search(center, cutoff)  # 5.0 for distance in angstrom
        residue_list = PDB.Selection.unfold_entities(neighbors, 'R')  # R for residues
        if len(residue_list) > 0:
            return True
    return False

def is_contact_CM(res_1, other_atoms, cutoff):
    center = 0
    for atom in res_1:
        center = center + atom.get_coord()
    center = center / len(res_1)
    ns = PDB.NeighborSearch(other_atoms)
    neighbors = ns.search(center, cutoff)  # 5.0 for distance in angstrom
    residue_list = PDB.Selection.unfold_entities(neighbors, 'R')  # R for residues
    if len(residue_list) > 0:
        return True

    return False

def get_contacts(residue_map, all_atoms, verbose, cutoff):
	progress = 0
	contacts = []
	for residue in residue_map:
		progress+=1
		if len(verbose)>0:
			print verbose,progress,"out of",len(residue_map)
		atom_list = residue_map[residue]

        # select a subset of all_atoms that is distance between atom_list and all_atoms are less than 4.5
		outcome = is_contact(atom_list,all_atoms,cutoff)
		if outcome:
			contacts.append(residue)
	return contacts


#Filter out all the atoms from the chain,residue map given by residue_map
def get_all_atoms(residue_map):
	all_atoms_out = []
	for residue in residue_map:
		for atom in residue_map[residue]:
			all_atoms_out.append(atom)
			#Set the b-factor to zero for coloring by contacts
			atom.set_bfactor(0.0)
	return all_atoms_out

# Get interface and contact residues without writing into a file
def get_target_residues(interface, contacts):
    Res = []
    for elem in interface:
        splitted = elem.split("_")
        resname = str((splitted[1]+splitted[2]).strip())
        chain = splitted[0]
        #contact or neighbor of interface?
        coninf = "I" #Interface neighbor
        if elem in contacts:
            coninf = "C" #Contact
        Res.append([chain, resname, coninf])

    return Res


#----------------------------------------------------------
# Get sequence information a PDB structure with Given Chain
#----------------------------------------------------------
def get_chSeq_df(pdb_file=None, pdb_ch=None):
    #----------------------------------------------------------
    Res3, Res1 = pdb_read_seq(pdb_file)
    aaSeq = []
    aaIdx = []
    aaCh = []
    for i_ch in pdb_ch:
        seqTmp = [aa[0] for aa in Res1[i_ch]]
        idxTmp = [aa[1] for aa in Res1[i_ch]]
        chTmp = len(seqTmp) * [i_ch]

        aaSeq = aaSeq + seqTmp
        aaIdx = aaIdx + idxTmp
        aaCh  = aaCh + chTmp

    seqInfo = pd.DataFrame({'SEQ':aaSeq, 'INDEX':aaIdx, 'CHAIN':aaCh})

    return seqInfo

def get_contact_labels(pdb_file1=None, pdb_file2=None, pdb_ch1=None, pdb_ch2=None, cutoff=4.5):
    seqInfo1 = get_chSeq_df(pdb_file=pdb_file1, pdb_ch=pdb_ch1)
    seqInfo2 = get_chSeq_df(pdb_file=pdb_file2, pdb_ch=pdb_ch2)

    #-- read PDB
    pdbParser = PDB.PDBParser()
    s1 = pdbParser.get_structure('structure1', pdb_file1)
    s2 = pdbParser.get_structure('structure2', pdb_file2)

    #-- get hierarhical structure
    # 'A'=Atom, 'R'=Residue, 'C' =Chain, 'M'=Model, 'S'=Structure
    # Chains1[0].id or tmp[0]_id will return the name of Chain
    Chains1 = PDB.Selection.unfold_entities(s1, 'C') # C for chains
    Chains2 = PDB.Selection.unfold_entities(s2, 'C') # C for chains

    #get the mapping from chain,residue id to the atom lists
    # atoms_1 contaions all atoms in a structure
    # input_1 and input_2 build 'residue_map'
    ResMap1 = get_atom_list(Chains1, pdb_ch1)  # choose the chain to be searched for contacts (e.g., chain 'A')
    ResMap2 = get_atom_list(Chains2, pdb_ch2)  # choose the chain to be searched for contacts (e.g., chain 'A')

    #get all Atoms of each residue map
    Atoms1 = get_all_atoms(ResMap1)
    Atoms2 = get_all_atoms(ResMap2)

    # Get contact residues
    #run neighbor search on both instances - not optimal but good enough for most imaginable applications.
    # ** ResMap1 against all atoms in S2 (i.e, Atoms2)
    conRes1 = get_contacts(ResMap1, Atoms2, '', cutoff)
    conRes2 = get_contacts(ResMap2, Atoms1, '', cutoff)
    #conRes1 = get_contacts(ResMap1, Atoms2, "First PDB, residue ", cutoff)
    #conRes2 = get_contacts(ResMap2, Atoms1, "Second PDB, residue ", cutoff)

    conRes1 = pd.DataFrame([res.split("_") for res in conRes1], columns=['CHAIN','INDEX','ICODE'])
    conRes2 = pd.DataFrame([res.split("_") for res in conRes2], columns=['CHAIN','INDEX','ICODE'])

    # get_chSeq_df() returns 'int' type and conRes get_contacts() returns 'str(object)' for 'INDEX' column,
    # so 'inner' cannot be applied between different data types. Convert conRes1 to 'int' type
    # conRes1['INDEX'] = conRes1['INDEX'].apply(pd.to_numeric)
    # conRes2['INDEX'] = conRes2['INDEX'].apply(pd.to_numeric)
    # mtxData = pd.merge(seqInfo1, conRes1, how='inner', on=['CHAIN', 'INDEX'])
    CLASS = [0] * seqInfo1.shape[0]
    for i in range(conRes1.shape[0]):
        i_ch = conRes1['CHAIN'][i]
        i_idx = conRes1['INDEX'][i]
        tg_idx = seqInfo1.index[(seqInfo1['CHAIN']==i_ch) & (seqInfo1['INDEX']== int(i_idx))].tolist()
        if len(tg_idx) > 0:
            CLASS[tg_idx[0]] = 1

    seqInfo1['CLASS'] = CLASS

    CLASS = [0] * seqInfo2.shape[0]
    for i in range(conRes2.shape[0]):
        i_ch = conRes2['CHAIN'][i]
        i_idx = conRes2['INDEX'][i]
        tg_idx = seqInfo2.index[(seqInfo2['CHAIN'] == i_ch) & (seqInfo2['INDEX'] == int(i_idx))].tolist()
        if len(tg_idx) > 0:
            CLASS[tg_idx[0]] = 1

    seqInfo2['CLASS'] = CLASS

    return seqInfo1, seqInfo2

def get_all_labels(pdb_file1=None, pdb_file2=None, pdb_ch1=None, pdb_ch2=None, cutoff=4.5, i_cutoff=5.0):
    seqInfo1 = get_chSeq_df(pdb_file=pdb_file1, pdb_ch=pdb_ch1)
    seqInfo2 = get_chSeq_df(pdb_file=pdb_file2, pdb_ch=pdb_ch2)

    #-- read PDB
    pdbParser = PDB.PDBParser()
    s1 = pdbParser.get_structure('structure1', pdb_file1)
    s2 = pdbParser.get_structure('structure2', pdb_file2)

    #-- get hierarhical structure
    # 'A'=Atom, 'R'=Residue, 'C' =Chain, 'M'=Model, 'S'=Structure
    # Chains1[0].id or tmp[0]_id will return the name of Chain
    Chains1 = PDB.Selection.unfold_entities(s1, 'C') # C for chains
    Chains2 = PDB.Selection.unfold_entities(s2, 'C') # C for chains

    #get the mapping from chain,residue id to the atom lists
    # atoms_1 contaions all atoms in a structure
    # input_1 and input_2 build 'residue_map'
    ResMap1 = get_atom_list(Chains1, pdb_ch1)  # choose the chain to be searched for contacts (e.g., chain 'A')
    ResMap2 = get_atom_list(Chains2, pdb_ch2)  # choose the chain to be searched for contacts (e.g., chain 'A')

    #get all Atoms of each residue map
    Atoms1 = get_all_atoms(ResMap1)
    Atoms2 = get_all_atoms(ResMap2)

    # Get contact residues
    #run neighbor search on both instances - not optimal but good enough for most imaginable applications.
    # ** ResMap1 against all atoms in S2 (i.e, Atoms2)
    conRes1 = get_contacts(ResMap1, Atoms2, '', cutoff)
    conRes2 = get_contacts(ResMap2, Atoms1, '', cutoff)
    #conRes1 = get_contacts(ResMap1, Atoms2, "First PDB, residue ", cutoff)
    #conRes2 = get_contacts(ResMap2, Atoms1, "Second PDB, residue ", cutoff)

    DF_conRes1 = pd.DataFrame([res.split("_") for res in conRes1], columns=['CHAIN','INDEX','ICODE'])
    DF_conRes2 = pd.DataFrame([res.split("_") for res in conRes2], columns=['CHAIN','INDEX','ICODE'])

    # get_chSeq_df() returns 'int' type and conRes get_contacts() returns 'str(object)' for 'INDEX' column,
    # so 'inner' cannot be applied between different data types. Convert conRes1 to 'int' type
    # conRes1['INDEX'] = conRes1['INDEX'].apply(pd.to_numeric)
    # conRes2['INDEX'] = conRes2['INDEX'].apply(pd.to_numeric)
    # mtxData = pd.merge(seqInfo1, conRes1, how='inner', on=['CHAIN', 'INDEX'])
    CONTACT = [0] * seqInfo1.shape[0]
    for i in range(DF_conRes1.shape[0]):
        i_ch = DF_conRes1['CHAIN'][i]
        i_idx = DF_conRes1['INDEX'][i]
        tg_idx = seqInfo1.index[(seqInfo1['CHAIN'] == i_ch) & (seqInfo1['INDEX'] == int(i_idx))].tolist()
        if len(tg_idx) > 0:
            CONTACT[tg_idx[0]] = 1

    seqInfo1['CONTACT'] = CONTACT

    CONTACT = [0] * seqInfo2.shape[0]
    for i in range(DF_conRes2.shape[0]):
        i_ch = DF_conRes2['CHAIN'][i]
        i_idx = DF_conRes2['INDEX'][i]
        tg_idx = seqInfo2.index[(seqInfo2['CHAIN'] == i_ch) & (seqInfo2['INDEX'] == int(i_idx))].tolist()
        if len(tg_idx) > 0:
            CONTACT[tg_idx[0]] = 1

    seqInfo2['CONTACT'] = CONTACT

    # Get interface residues
    # run neighbor search on both instances - not optimal but good enough for most imaginable applications.
    # ** ResMap1 against all atoms in S1 (i.e, Atoms1)
    contact_map_1 = []
    for residue in conRes1:
        for atom in ResMap1[residue]:
            atom.set_bfactor(100.0)
            contact_map_1.append(atom)

    contact_map_2 = []
    for residue in conRes2:
        for atom in ResMap2[residue]:
            atom.set_bfactor(100.0)
            contact_map_2.append(atom)

    intRes1 = conRes1
    intRes2 = conRes2
    if i_cutoff > 0:
        intRes1 = get_contacts(ResMap1, contact_map_1, '', i_cutoff)
        intRes2 = get_contacts(ResMap2, contact_map_2, '', i_cutoff)
        #intRes1 = get_contacts(ResMap1, contact_map_1, "First PDB, interfacial residue ", i_cutoff)
        #intRes2 = get_contacts(ResMap2, contact_map_2, "Second PDB, interfacial residue ", i_cutoff)

    DF_intRes1 = pd.DataFrame([res.split("_") for res in intRes1], columns=['CHAIN', 'INDEX', 'ICODE'])
    DF_intRes2 = pd.DataFrame([res.split("_") for res in intRes2], columns=['CHAIN', 'INDEX', 'ICODE'])

    INTERFACE = [0] * seqInfo1.shape[0]
    for i in range(DF_intRes1.shape[0]):
        i_ch = DF_intRes1['CHAIN'][i]
        i_idx = DF_intRes1['INDEX'][i]
        tg_idx = seqInfo1.index[(seqInfo1['CHAIN'] == i_ch) & (seqInfo1['INDEX'] == int(i_idx))].tolist()
        if len(tg_idx) > 0:
            INTERFACE[tg_idx[0]] = 1
    seqInfo1['INTERFACE'] = INTERFACE

    INTERFACE = [0] * seqInfo2.shape[0]
    for i in range(DF_intRes2.shape[0]):
        i_ch = DF_intRes2['CHAIN'][i]
        i_idx = DF_intRes2['INDEX'][i]
        tg_idx = seqInfo2.index[(seqInfo2['CHAIN'] == i_ch) & (seqInfo2['INDEX'] == int(i_idx))].tolist()
        if len(tg_idx) > 0:
            INTERFACE[tg_idx[0]] = 1
    seqInfo2['INTERFACE'] = INTERFACE

    return seqInfo1, seqInfo2

#----------------------------------------------------------
# Download PDB files and store them into local HDD with PDB format
#----------------------------------------------------------
def downloadPDB(pdb_id, PDB_HOME):
    pdbDL = PDB.PDBList()
    cif_file = '{}/{}.cif'.format(PDB_HOME, pdb_id)
    while not os.path.exists(cif_file):
        try: 
            cif_file = pdbDL.retrieve_pdb_file(pdb_code=pdb_id, pdir=PDB_HOME, file_format='mmCif',  obsolete=False)
        except:
            pass
    
    # read any file format
    #cif_file = "{}/{}.cif".format(PDB_HOME, pdb_id.lower())
    cifParser = PDB.MMCIFParser()
    structure = cifParser.get_structure(pdb_id, cif_file)

    pdb_file = "{}/{}.pdb".format(PDB_HOME, pdb_id.lower())
    pdbIO = PDB.PDBIO()
    pdbIO.set_structure(structure[0])
    pdbIO.save(pdb_file)
    return 0


#----------------------------------------------------------
# Get all possible combinations of pairs of chains for 3 Chains
#----------------------------------------------------------
def getPairsL3(vec):
    Pool = []
    #-- unique two pairs 1:1
    tmp = list(itertools.combinations(vec, 2))
    for i in tmp:
        Pool.append(list(i))
    
        #-- unique 1:2
    for i in range(len(vec)):
        newV = copy.copy(vec)
        newV.remove(vec[i])
        tmp = list(itertools.combinations(newV, 2))
        for j in tmp:
            Chains = ''.join(j)
            tmpPair = [vec[i], Chains]
            tmpPair.sort()
            Pool.append(tmpPair)
    
    #-- make unique list of lists and sorting them
    outData = [list(x) for x in set(tuple(x) for x in Pool)]
    outData = sorted(outData, key=operator.itemgetter(0))
    return outData



#----------------------------------------------------------
# Get all possible combinations of pairs of chains for 4 Chains
#----------------------------------------------------------
def getPairsL4(vec):
    Pool = []
    #-- unique two pairs 1:1
    tmp = list(itertools.combinations(vec, 2))
    for i in tmp:
        Pool.append(list(i))
    
        #-- unique 1:2
    for i in range(len(vec)):
        newV = copy.copy(vec)
        newV.remove(vec[i])
        tmp = list(itertools.combinations(newV, 2))
        for j in tmp:
            Chains = ''.join(j)
            tmpPair = [vec[i], Chains]
            tmpPair.sort()
            Pool.append(tmpPair)
            
    #-- unique 2:2
    tmp = list(itertools.combinations(vec, 2))
    for i in tmp:
        tmp2 = filter(lambda x: not(x in list(i)), vec)
        tmpPair = [''.join(list(i)), ''.join(tmp2)]
        tmpPair.sort()
        Pool.append(tmpPair)
                     
    #-- unique 1:3
    for i in range(len(vec)):
        newV = copy.copy(vec)
        newV.remove(vec[i])
        tmp = list(itertools.combinations(newV, 3))
        for j in tmp:
            Chains = ''.join(j)
            tmpPair = [vec[i], Chains]
            tmpPair.sort()
            Pool.append(tmpPair)
    
    #-- make unique list of lists and sorting them
    outData = [list(x) for x in set(tuple(x) for x in Pool)]
    outData = sorted(outData, key=operator.itemgetter(0))
    return outData


#----------------------------------------------------------
# Count Motif
# E.g.,
#   Seq = 'ANLDRSNDKVYENVTGLV'
#   Motif = ['DR', 'ND', 'VT']
# ** Motif's length must be SAME: ['DR', 'ND']->(O) ['DRT', 'ND']->(X)
#----------------------------------------------------------
def countMotif(Seq, Motif):
    win_size = len(Motif[0])
    obj = jc_window.winSeq(Seq)
    objWin = obj.makelinear(win_size=win_size, intv=2)
    winSeq = jc_window.retrieveSeq(objWin, len(Seq))
    sortWin = map(lambda x: ''.join(sorted(list(x))), winSeq)
    cntMotif = collections.Counter(sortWin)

    cntVec = [0] * len(Motif)
    Keys = cntMotif.keys()
    for k in Keys:
        if k in Motif:
            cntVec[Motif.index(k)] = cntMotif[k]

    return Motif, cntVec

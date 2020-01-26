import jctools
import pprint
import pandas as pd
from Bio import SeqIO
from Bio import PDB

from Bio.PDB import PDBIO
import os
import random
import ntpath

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


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


#----------------------------------------------------------
# Read Chains from PDB
#----------------------------------------------------------
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
# Get interface and contact lists
#----------------------------------------------------------
pdb_file1 = '/Users/jjeong/local/project_dev/ppi/codes/utils/1A2Y.pdb'
pdb_file2 = '/Users/jjeong/local/project_dev/ppi/codes/utils/1A2Y.pdb'
pdb_id1 = '1A2Y'
pdb_id2 = '1A2Y'
pdb_ch1 = 'AB'
pdb_ch2 = 'C'


#----------------------------------------------------------
# Get sequence information a PDB structure with Given Chain
#----------------------------------------------------------
def get_chSeq_df(pdb_file=None, pdb_ch=None):
    #----------------------------------------------------------
    Res3, Res1 = jctools.pdb_read_seq(pdb_file)
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
    s1 = pdbParser.get_structure(pdb_id1, pdb_file1)
    s2 = pdbParser.get_structure(pdb_id2, pdb_file2)

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
    s1 = pdbParser.get_structure(pdb_id1, pdb_file1)
    s2 = pdbParser.get_structure(pdb_id2, pdb_file2)

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


#cS1, cS2 = get_contact_labels(pdb_file1=pdb_file1, pdb_file2=pdb_file2, pdb_ch1=pdb_ch1, pdb_ch2=pdb_ch2, cutoff=4.5)
S1, S2 = jctools.get_all_labels(pdb_file1=pdb_file1, pdb_file2=pdb_file2, pdb_ch1=pdb_ch1, pdb_ch2=pdb_ch2, cutoff=4.5, i_cutoff=5.0)

"""
# Get interface residues
#run neighbor search on both instances - not optimal but good enough for most imaginable applications.
# ** ResMap1 against all atoms in S1 (i.e, Atoms1)
i_cutoff=5.0
intRes1 = conRes1
intRes2 = conRes2
if i_cutoff > 0:
    intRes1 = get_contacts(ResMap1, Atoms1, "First PDB, interfacial residue ", i_cutoff)
    intRes1 = pd.DataFrame([res.split("_") for res in intRes1], columns=['CHAIN', 'INDEX', 'ICODE'])
    intRes2 = get_contacts(ResMap2, Atoms2, "Second PDB, interfacial residue ", i_cutoff)
    intRes2 = pd.DataFrame([res.split("_") for res in intRes2], columns=['CHAIN', 'INDEX', 'ICODE'])

"""




"""
# 'A'=Atom, 'R'=Residue, 'C' =Chain, 'M'=Model, 'S'=Structure
# tmp[0].id or tmp[0]_id will return values
Structure2 = PDB.Selection.unfold_entities(structure, 'C')
Chains2 = [c._id for c in Structure2]

tmp = PDB.Selection.unfold_entities(structure, 'R') # this will include none canonical 20 AAs e.g., 'HOH'
Res = [res.get_resname() for res in tmp]
Atoms = [atom.get_list() for atom in tmp]

Chains = pdb_read_chain(pdb_file)
residue_map = get_atom_list(Structure2, Chains2)

all_atoms = get_all_atoms(residue_map)
contacts_1 = get_contacts(residue_map, all_atoms, "First molecule, residue ", 4.5)
"""
print(pdb_id)
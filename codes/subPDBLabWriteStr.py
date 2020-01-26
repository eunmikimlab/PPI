"""
   ::Purpose::
   Dry run to extract PDB Chains corresponding to PDB name and write them into a file
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
import itertools
import copy
import operator

#------------------------------------------------------------------------------
# Parameters
# python PDBChains.py -i 14gs
# -o /home/jj140/scratch/ppi/batch/EC2/outputs
# -r /home/jj140/scratch/ppi/dbs/pdb_chain_enzyme.csv
# -s 10078
# -e 40899
#------------------------------------------------------------------------------
python = '/home/jj140/local/bin/anaconda2/envs/ppi/bin/python'          # Python executable 
code_path = '/home/jj140/project/ppi/codes/PDBLabWriteStr.py'              # Path of a code to be run
qName = 'short'                                                         # Name of partition (queue) to be submitted 
constraint='scratch2'                                                   # set to None if N/A - running codes under particular File system
cont_dist = 4.5                                                         # distant to decide contact residues - direct contact to another protein
intf_dist = 0                                                        # distant to decide interface residue - residues that are near (e.g., 10.0A) contact residues

# python PDBLabWriteStr.py --f1 /home/jj140/scratch/ppi/batch/pdb/pdbOnly/10gs.pdb
# --f2 /home/jj140/scratch/ppi/batch/pdb/pdbOnly/10gs.pdb
# --c1 A --c2 B --c 4.5 --i 10.0
# --jobid /home/jj140/scratch/ppi/batch/PDBlabels
dbName = 'SCOP'
HOME_DIR = '/home/jj140/scratch/ppi/batch/PDBlabels/{}'.format(dbName)      # Home directory where data want to be stored
pdbFile = '/home/jj140/scratch/ppi/batch/{}_mergedDBs.csv'.format(dbName)   # PDB and all associated chains
pdbDir = '/home/jj140/scratch/ppi/batch/pdb/pdbOnly'                        # source of PDB files 

s_idx = int(0)
n_job = int(10)


#------------------------------------------------------------------------------
# Define path where data are stored
#------------------------------------------------------------------------------
OUT_DIR = '{}'.format(HOME_DIR)
try:
    os.makedirs(OUT_DIR)
except OSError:
    pass

out_script = '{}/scripts'.format(OUT_DIR)
try:
    os.makedirs(out_script)
except OSError:
    pass

out_result = '{}/outputs'.format(OUT_DIR)
try:
    os.makedirs(out_result)
except OSError:
    pass

out_log = '{}/slurm_log'.format(OUT_DIR)
try:
    os.makedirs(out_log)
except OSError:
    pass




#- get start and end index
def getSteps(s_idx, n_job, len_pdb):
    IDX = []
    e_tmp = 0
    while e_tmp <= len_pdb:
        if e_tmp == 0:
            s_tmp = s_idx
            e_tmp = 1
        else:
            e_tmp = s_tmp + n_job
            IDX.append([s_tmp, e_tmp])
            s_tmp = e_tmp
    return IDX


#------------------------------------------------------------------------------
# Search PDB and extract corresponding Chains
#------------------------------------------------------------------------------
dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})
len_pdb = len(dfPDB)
Steps = getSteps(s_idx, n_job, len_pdb)

obj = jctools.jobSLURM()
scrFiles = []
for step in Steps:
    s_idx = step[0]
    e_idx = step[1]
    
    src_file = '{}/job{}_{}.sh'.format(out_script, s_idx, e_idx)
    scrFiles.append(src_file)
    
    selPDB = dfPDB['PDB'][s_idx:e_idx].tolist()
    jName = 'job{}'.format(s_idx)
    
    exe_code = []
    for pdb_id in selPDB:
        Chains = list(dfPDB.loc[dfPDB['PDB'] == pdb_id][dbName].tolist()[0])
        num_ch = len(Chains)
        
        Pools = []
        if num_ch == 2:
            Pools.append(Chains)
        elif num_ch == 3:
            Pools = jctools.getPairsL3(Chains)
        elif num_ch == 4:
            Pools = jctools.getPairsL4(Chains)
        else:
            Chains = Chains[0:4]
            Pools = jctools.getPairsL4(Chains)
        
        pdb_file_1 = '{}/{}.pdb'.format(pdbDir, pdb_id)
        pdb_file_2 = '{}/{}.pdb'.format(pdbDir, pdb_id)
        
        for pair in Pools:
            pdb_ch1 = pair[0]
            pdb_ch2 = pair[1]
            # python PDBLabWriteStr.py --f1 /home/jj140/scratch/ppi/batch/pdb/pdbOnly/10gs.pdb --f2 /home/jj140/scratch/ppi/batch/pdb/pdbOnly/10gs.pdb --c1 A --c2 B --c 4.5 --i 10.0 --jobid /home/jj140/scratch/ppi/batch/PDBlabels    
            cmd = '{} {} --f1 {} --f2 {} --c1 {} --c2 {} --c {} --i {} --jobid {}'.format(python, code_path, pdb_file_1, pdb_file_1, pdb_ch1, pdb_ch2, cont_dist, intf_dist, out_result)
            exe_code.append(cmd)
        
    strExec = '{}'.format('\n'.join(exe_code))
    obj.writejob(outFile=src_file, jName=jName, qName=qName, constraint=constraint, home_dir=HOME_DIR, exe_code=strExec, log_dir=out_log)
    

# Submit Jobs
job_summary = '{}/job_summary.txt'.format(OUT_DIR)
obj.submitjob(fileNames=scrFiles, home_dir=out_log, outFile=job_summary)

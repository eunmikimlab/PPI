"""
   ::Purpose::
   Write lists of dssp and hssp download list to use rsync
   rsync -auvh --progress --files-from=/path/to/files.txt / /destination/path/
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
#------------------------------------------------------------------------------
python = '/home/jj140/local/bin/anaconda2/envs/ppi/bin/python'          # Python executable 
qName = 'short'                                                         # Name of partition (queue) to be submitted 
constraint='scratch2'                                                   # set to None if N/A - running codes under particular File system

dbName = 'hssp'
prefix = 'rsync://rsync.cmbi.ru.nl'
#prefix = '/rsync.cmbi.ru.nl'
pdbFile = '/home/jj140/scratch/ppi/batch/EC_mergedDBs.csv'.format(dbName)   # PDB and all associated chains
HOME_DIR = '/home/jj140/scratch/ppi/batch/{}'.format(dbName)      # Home directory where data want to be stored

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
        if dbName == 'dssp':
                cmd = '\nRC=1\n'
                cmd = cmd + 'while [[ $RC -ne 0 ]]\n'
                cmd = cmd + 'do\n'
                cmd = cmd + '\t' + 'rsync -auvh {}/dssp/{}.dssp {}/.\n'.format(prefix, pdb_id.lower(), out_result)
                cmd = cmd + '\t' + 'RC=$?\n'
                cmd = cmd + 'done\n'
                
        elif dbName == 'hssp':
                cmd = '\nRC=1\n'
                cmd = cmd + 'while [[ $RC -ne 0 ]]\n'
                cmd = cmd + 'do\n'
                cmd = cmd + '\t' + 'rsync -auvh {}/hssp/{}.hssp.bz2 {}/.\n'.format(prefix, pdb_id.lower(), out_result)
                cmd = cmd + '\t' + 'RC=$?\n'
                cmd = cmd + 'done\n'
                cmd = cmd + '\n' + 'bzip2 -d {}/{}.hssp.bz2'.format(out_result, pdb_id)

        exe_code.append(cmd)
        
    strExec = '{}'.format('\n'.join(exe_code))
    obj.writejob(outFile=src_file, jName=jName, qName=qName, constraint=constraint, home_dir=HOME_DIR, exe_code=strExec, log_dir=out_log)
    

# Submit Jobs
job_summary = '{}/job_summary.txt'.format(OUT_DIR)
obj.submitjob(fileNames=scrFiles, home_dir=out_log, outFile=job_summary)

"""
   ::Purpose::
   This code submit jobs to HPC to get PDB and Chain pairs
   This code is optimized for O2 Slurm cluster
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 13, 2018"

import jctools
import pandas as pd
import os

#------------------------------------------------------------------------------
# Parameters
#------------------------------------------------------------------------------
#[O2]

python = '/home/jj140/local/bin/anaconda2/envs/ppi/bin/python'          # Python executable 
code_path = '/home/jj140/project/ppi/codes/PDBDownload.py'              # Path of a code to be run
HOME_DIR = '/home/jj140/scratch/ppi/batch/pdb'                          # Home directory where data want to be stored
qName = 'short'                                                         # Name of partition (queue) to be submitted 
constraint='scratch2'                                                   # set to None if N/A - running codes under particular File system
pdbFile = '/home/jj140/scratch/ppi/dbs/outputs/pdb_missing.csv'         # Unique PDB list - nx1 Column Head = 'PDB'

s_idx = int(0)
n_job = int(100)

#[OSX]
"""
python = '/Users/jjeong/local/programs/anaconda2/envs/ppi/bin/python'
code_path = '/Users/jjeong/local/project_dev/ppi/codes/PDBChains.py'
HOME_DIR = '/Users/jjeong/local/project_dev/ppi/slurm'
qName = 'short'
constraint='scratch2' # set to None if N/A
pdbFile = '/Users/jjeong/local/project_dev/ppi/dbs/outputs/unique_human_pdb.csv'
dbFile  = '/Users/jjeong/local/project_dev/ppi/dbs/pdb_chain_enzyme.csv'
"""


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
# Main codes
#------------------------------------------------------------------------------
dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})
len_pdb = len(dfPDB)

#s_idx = 3
#n_job = 4
#len_pdb = 23

Steps = getSteps(s_idx, n_job, len_pdb)


obj = jctools.jobSLURM()
scrFiles = []
for step in Steps:
    s_idx = step[0]
    e_idx = step[1]

    src_file = '{}/job{}_{}.sh'.format(out_script, s_idx, e_idx)

    scrFiles.append(src_file)

    jName = 'job{}'.format(s_idx)
    exe_code = '{} {} -i {} -o {} -s {} -e {}'.format(python, code_path, pdbFile, out_result, s_idx, e_idx)
    obj.writejob(outFile=src_file, jName=jName, qName=qName, constraint=constraint, home_dir=HOME_DIR, exe_code=exe_code, log_dir=out_log)
    

# Submit Jobs
job_summary = '{}/job_summary.txt'.format(OUT_DIR)
obj.submitjob(fileNames=scrFiles, home_dir=out_log, outFile=job_summary)




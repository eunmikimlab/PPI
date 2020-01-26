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
"""
python = '/home/jj140/local/bin/anaconda2/envs/ppi/bin/python'
code_path = '/home/jj140/project/ppi/codes/PDBChains.py'
HOME_DIR = '/home/jj140/scratch/ppi/batch'
qName = 'short'
constraint='scratch2' # set to None if N/A
pdbFile = '/Users/jjeong/local/project_dev/ppi/dbs/outputs/unique_human_pdb.csv'
dbFile  = '/Users/jjeong/local/project_dev/ppi/dbs/pdb_chain_enzyme.csv'
"""

#[OSX]
python = '/Users/jjeong/local/programs/anaconda2/envs/ppi/bin/python'
code_path = '/Users/jjeong/local/project_dev/ppi/codes/PDBChains.py'
HOME_DIR = '/Users/jjeong/local/project_dev/ppi/slurm'
qName = 'short'
constraint='scratch2' # set to None if N/A
pdbFile = '/Users/jjeong/local/project_dev/ppi/dbs/outputs/unique_human_pdb.csv'
dbFile  = '/Users/jjeong/local/project_dev/ppi/dbs/pdb_chain_enzyme.csv'

OUT_DIR = '{}/pdb_chains'.format(HOME_DIR)
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

#------------------------------------------------------------------------------
# Main codes
#------------------------------------------------------------------------------

dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})
#len_pdb = len(dfPDB)
len_pdb = 3
obj = jctools.jobSLURM()
scrFiles = []
for i in range(len_pdb):
    pdb_id = dfPDB['PDB'][i]
    src_file = '{}/job{}_{}.sh'.format(out_script, i, pdb_id)
    out_file = '{}/{}_{}.txt'.format(out_result, i, pdb_id)

    scrFiles.append(src_file)

    jName = '{}'.format(pdb_id)
    exe_code = '{} {} -i 14gs -o {} -r {}'.format(python, code_path, out_file, dbFile)
    obj.writejob(outFile=src_file, jName=jName, qName=qName, constraint=constraint, home_dir=HOME_DIR, exe_code=exe_code, log_dir=out_log)

obj.submitjob(fileNames=scrFiles, home_dir=out_log)

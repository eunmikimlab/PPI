"""
   ::Purpose::
   Extract PDB Chains corresponding to PDB name and write them into a file
   This code is optimized for O2 Slurm cluster
"""
__author__      = "Jong Cheol Jeong"
__copyright__   = "Copyright 2018, Jeong Lab."
__license__     = "Jeong Lab."
__email__       = "JongCheol.Jeong@uky.edu"
__date__        = "Sep 13, 2018"


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

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help='A file containing the list of PDB IDs and associated Chains', action='store', required=True)
parser.add_argument("-o", "--outdir", help='Output directory to store PDB and Chain pairs', action='store', required=True)
parser.add_argument("-d", "--dbname", help='Name of DB used as column name', action='store', required=True)
parser.add_argument("-s", "--start", help='The index of PDB data to be started', action='store', required=True)
parser.add_argument("-e", "--end", help='The index of PDB data to be ended', action='store', required=True)
parser.add_argument("-p", "--pdbhome", help='PDB home directory', action='store', default='/home/jj140/scratch/ppi/batch/pdb/pdbOnly')

args = parser.parse_args()

pdbFile = args.input
outDir = args.outdir
dbName = args.dbname
pdbDir = args.pdbhome

s_idx = int(args.start)
e_idx = int(args.end)


report_dir = '{}/report'.format(outDir)
try:
    os.makedirs(report_dir)
except OSError:
    pass

"""
#python PDBLabels.py -i /home/jj140/scratch/ppi/batch/mergedDBs.csv -o /home/jj140/scratch/ppi/batch/labels -d SCOP -s 10 -e 11
pdbFile = '/home/jj140/scratch/ppi/batch/mergedDBs.csv'
outDir = '/home/jj140/scratch/ppi/batch/labels'
dbName = 'EC'
pdbDir = '/home/jj140/scratch/ppi/batch/pdb/pdbOnly'
s_idx = 0
e_idx = 1
"""
#------------------------------------------------------------------------------
# Search PDB and extract corresponding Chains
#------------------------------------------------------------------------------
dfPDB = pd.read_csv(pdbFile, comment='#', dtype={'PDB': object})
selPDB = dfPDB['PDB'][s_idx:e_idx].tolist()

out_report = '{}/report_{}.tsv'.format(report_dir, s_idx)
with open(out_report, 'w') as f:
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
            
            outFile1 = '{}/{}_{}-{}_{}.csv'.format(outDir, pdb_id, pdb_ch1, pdb_ch2, pdb_ch1)
            outFile2 = '{}/{}_{}-{}_{}.csv'.format(outDir, pdb_id, pdb_ch1, pdb_ch2, pdb_ch2)
            #S1, S2 = jctools.get_all_labels(pdb_file1=pdb_file_1, pdb_file2=pdb_file_2, pdb_ch1=pdb_ch1, pdb_ch2=pdb_ch2, cutoff=4.5, i_cutoff=5.0)
            S1, S2 = jctools.get_contact_labels(pdb_file1=pdb_file_1, pdb_file2=pdb_file_2, pdb_ch1=pdb_ch1, pdb_ch2=pdb_ch2, cutoff=4.5)
            S1.to_csv(outFile1, index=False)
            S2.to_csv(outFile2, index=False)
            
            strReport = '{}\t{}\t{}\t{}\t{}\n'.format(pdb_id, pdb_ch1, pdb_ch2, outFile1, outFile2)
            f.write(strReport)


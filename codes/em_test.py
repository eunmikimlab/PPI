import jctools
import pandas as pd

pdb_file1 = '/Users/emkim/local/projects/ppi/pdb/1fvc.pdb'
pdb_file2 = '/Users/emkim/local/projects/ppi/pdb/1fvc.pdb'
pdb_ch1 = 'AB'
pdb_ch2 = 'CD'
Chains = jctools.pdb_read_chain(pdb_file1)

S1, S2 = jctools.get_all_labels(pdb_file1=pdb_file1, pdb_file2=pdb_file2, pdb_ch1=pdb_ch1, pdb_ch2=pdb_ch2, cutoff=4.5, i_cutoff=5.0)

S1.to_csv('/Users/emkim/local/projects/ppi/output/1fvc_AB.csv')
S2.to_csv('/Users/emkim/local/projects/ppi/output/1fvc_CD.csv')

print('End!')
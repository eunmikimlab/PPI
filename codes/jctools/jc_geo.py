"""
   Jong Cheol's Package contains GEO related jobs
"""
import re

class geo:
    def __init__(self):
	self.home_dir = '~';
    
    '''-----------------------------
	Parsing the TXT file downloaded from https://www.ncbi.nlm.nih.gov/gds
	It contains following information:
	    <example>
	    Organism:	Caenorhabditis elegans
	    Source name:	C. elegans
	    Platform: GPL19757
	    Series: GSE81707 GSE81708 
	    FTP download: 
	    SRA Run Selector: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRX1787486
	    Accession: GSM2171673
	    ID: 302171673
       -----------------------------'''
       
    def findIdx(self, data, y):
	outData = [i for i,x in enumerate(data) if x == y];
	return outData
	
    def gds_text(self, inFile=None):
	Templates = ['organism',
		     'source name',
		     'platform',
		     'series',
		     'ftp download',
		     'sra run selector',
		     'accession',
		     'id'];
	

	outList = [];
	
	with open(inFile, 'r') as f:
	    tmpVal  = [0 for i in range(len(Templates))];
	    for line in f:
		for i in range(len(Templates)):
		    ptrn = '.*{}:(.*)'.format(Templates[i])
		    rst = re.match(ptrn, line, re.IGNORECASE);
		    
		    if rst is not None:
			tmpVal[i] = rst.group(1).strip();
			
		    if len(self.findIdx(tmpVal, 0)) < 1:
			tmpVal[2] = re.split(' |\t', tmpVal[2])[0]
			tmpVal[6] = re.split(' |\t', tmpVal[6])[0]
			outList.append(tmpVal);
			#print("\n".join(tmpVal));
			#print '========================================='
			tmpVal  = [0 for i in range(len(Templates))];
			break;
	    f.closed
	    
	    return outList;
	    
	


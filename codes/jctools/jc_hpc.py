"""
   Jong Cheol's Package contain HPC related jobs
"""
import os
import subprocess
import time
import re


#------------------------------------------------------------------------------
# Job submission for SLURM
# Example Codes
# import jc_hpc
#
# obj = jc_hpc.jobSLURM()
# Files = []
#
# home_dir = '/home/jj140/scratch/ppi/batch'
# for i in range(10):
#   outname = './script/myfile{}.sh'.format(i)
#   Files.append(outname)
#   jName = 'JeongLab{}'.format(i)
#   qName = 'short'
#   constraint='scratch2'
#   exe_code = 'echo job_number={} > test{}.txt'.format(i,i)
#   obj.writejob(outFile=outname, jName=jName, qName=qName, constraint=constraint, home_dir=home_dir, exe_code=exe_code)
# bj.submitjob(fileNames=Files, home_dir=home_dir)
#------------------------------------------------------------------------------
class jobSLURM:
    def __init__(self):
        self.home_dir = '~'
    
    '''-----------------------------
    write bsub script
    -----------------------------'''
    def writejob(self, outFile=None,
           jName='JeongLab',
           qName='short',
           wTime='0-12:00',
           mem='4G',
           nCore='1',
           constraint=None,
           home_dir='~/.',
           exe_code=None,
           env_code='~/.bashrc',
           log_dir=None):

        cmt = '#!/bin/bash\n\n'
        cmt = '{}#SBATCH -J {}\n'.format(cmt, jName)
        cmt = '{}#SBATCH -p {}\n'.format(cmt, qName)
        cmt = '{}#SBATCH -t {}\n'.format(cmt, wTime)
        cmt = '{}#SBATCH -c {}\n'.format(cmt, nCore)

        if (constraint is not None):
            cmt = '{}#SBATCH --constraint="{}"\n'.format(cmt, constraint)

        if (mem is not None):
            cmt = '{}#SBATCH --mem={}\n'.format(cmt, mem)

        #-- build error directory
        if (log_dir is not None):
            slurm_log = log_dir
        else:
            slurm_log = '{}/slurm_log'.format(home_dir)

        try:
            os.makedirs(slurm_log)
        except OSError:
            pass

        cmt = '{}#SBATCH -o {}/bLab.%J.out\n'.format(cmt, slurm_log)
        cmt = '{}#SBATCH -e {}/bLab.%J.err\n\n'.format(cmt, slurm_log)

        cmt = '{}source {}\n'.format(cmt, env_code)
        cmt = '{}cd {}\n'.format(cmt, home_dir)
        cmt = '{}{}\n'.format(cmt, exe_code)

        if (outFile is not None):
            with open(outFile, 'w') as f:
                f.write(cmt)

        return cmt


    
    '''-----------------------------
    submit bsub script
    -----------------------------'''
    def submitjob(self,fileNames, home_dir='~/.', fps=2, elps_sec=0, outFile='./SLURM_jobLists.txt', verbose=True):
    
        #// creating 'bLab_lsf_err_log' directory to store '.err' and '.out' file
        err_dir = "{}/slurm_log".format(home_dir);
        out_dir = "{}/slurm_log".format(home_dir);

        if not os.path.exists(err_dir):
            os.makedirs(err_dir);

        if not os.path.exists(out_dir):
            os.makedirs(out_dir);

        cnt = 0;
        num_row = len(fileNames);
        outDATA = [[0 for c in range(2)] for r in range(num_row)];

        for i in fileNames:
            #// Submit jobs to LSF queue
            bsub = 'sbatch {}'.format(i);

            job_id = subprocess.check_output(bsub, shell=True)
            job_id = re.match(r"Submitted.+ (\d+)$", job_id)
            job_id = job_id.group(1);

            #// store submitted job information
            outDATA[cnt][0] = job_id;
            outDATA[cnt][1] = i

            if verbose:
                cmt = '\t[{}/{}] {}\t{}'.format(cnt+1, num_row, job_id, i)
                print(cmt)

            #// pause ELSP_SEC after sending jobs
            if (cnt % fps == 0):
                time.sleep(elps_sec)

            cnt += 1;

        with open(outFile, 'w') as f:
            f.write('JOB_ID\tFILE_NAME\n')
            for j in range(num_row):
                tmp = '{}\t{}\n'.format(outDATA[j][0], outDATA[j][1])
                f.write(tmp)
        f.closed



#------------------------------------------------------------------------------
# Job submission for LSF
#------------------------------------------------------------------------------
class jobLSF:
    def __init__(self):
        self.home_dir = '~';
    
    '''-----------------------------
    write bsub script
       -----------------------------'''
    def writejob(self, outFile=None,
           jName='bLab',
           qName='short',
           wTime='720',
           mem=None,
           tmp=None,
           nCore='1',
           priority='50',
           home_dir='~/.',
           exe_code='',
           env_code='~/.bashrc'):

        cmt = '#!/bin/bash\n\n';
        cmt = '{}#BSUB -J {}\n'.format(cmt, jName);
        cmt = '{}#BSUB -q {}\n'.format(cmt, qName);
        cmt = '{}#BSUB -W {}\n'.format(cmt, wTime);
        cmt = '{}#BSUB -n {}\n'.format(cmt, nCore);
        cmt = '{}#BSUB -sp {}\n'.format(cmt, priority);

        if (tmp is not None):
            cmt = '{}#BSUB -R "rusage[tmp={}]"\n'.format(cmt, tmp);

        if (mem is not None):
            cmt = '{}#BSUB -R "rusage[mem={}]"\n'.format(cmt, mem);

        cmt = '{}#BSUB -o ./bLab_lsf_out_log/bLab.%J.out\n'.format(cmt);
        cmt = '{}#BSUB -e ./bLab_lsf_err_log/bLab.%J.err\n\n'.format(cmt);

        cmt = '{}source {}\n'.format(cmt, env_code);
        cmt = '{}cd {}\n'.format(cmt, home_dir);
        cmt = '{}{}\n'.format(cmt, exe_code);

        if (outFile is not None):
            with open(outFile, 'w') as f:
                f.write(cmt)

        return(cmt)


    
    '''-----------------------------
    submit bsub script
       -----------------------------'''
    def submitjob(self,fileNames, fps=2, elps_sec=0, outFile='./lsfScriptList.txt', verbose=True):
    
        #// creating 'bLab_lsf_err_log' directory to store '.err' and '.out' file
        err_dir = './bLab_lsf_err_log';
        out_dir = './bLab_lsf_out_log';

        if not os.path.exists(err_dir):
            os.makedirs(err_dir);

        if not os.path.exists(out_dir):
            os.makedirs(out_dir);

        cnt = 0;
        num_row = len(fileNames);
        outDATA = [[0 for c in range(2)] for r in range(num_row)];

        for i in fileNames:
            #// Submit jobs to LSF queue
            bsub = 'bsub < {}'.format(i);

            job_id = subprocess.check_output(bsub, shell=True)
            job_id = re.match(r".*<(\d+)>.+", job_id)
            job_id = job_id.group(1);

            #// store submitted job information
            outDATA[cnt][0] = job_id;
            outDATA[cnt][1] = i;

            if verbose:
                cmt = '\t[{}/{}] {}\t{}'.format(cnt+1, num_row, job_id, i)
                print(cmt)

            #// pause ELSP_SEC after sending jobs
            if (cnt % fps == 0):
                time.sleep(elps_sec);

            cnt += 1

        with open(outFile, 'w') as f:
            f.write('JOB_ID\tFILE_NAME\n')
            for j in range(num_row):
                tmp = '{}\t{}\n'.format(outDATA[j][0], outDATA[j][1])
                f.write(tmp)

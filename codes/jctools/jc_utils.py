"""
   Jong Cheol's Package contain various utilities 
"""
import os
import subprocess
import time
import re


class jcDirectory:
    def __init__(self, home_dir):
        self.home_dir = home_dir
        self.results = []

    #-- get all subdirectories
    def getsubdir(self, Dirs=None, cnt=0):
        if Dirs is None:
            Dirs = os.walk(self.home_dir).next()[1]
            Results = map(lambda x: '{}/{}'.format(self.home_dir, x), Dirs)
            Dirs = Results
            self.results += Results
            self.getsubdir(Dirs, cnt)
            return self.results
        elif isinstance(Dirs, list):
            if len(Dirs) > 0:
                tmp_dir = os.walk(Dirs[0]).next()[1]
                if len(tmp_dir) > 0:
                    tmp_dir = map(lambda x: '{}/{}'.format(Dirs[0], x), tmp_dir)
                    Dirs += tmp_dir
                    self.results += tmp_dir
                cnt += 1
                del Dirs[0]
                #print '=========RESULTS========'
                for x in self.results: print x

                #print '=========DIRS========'
                for x in Dirs: print x
                self.getsubdir(Dirs, cnt)
                return self.results
            else:
                return self.results
        else:
            return self.results

    # -- get all files
    def listFiles(self, Dirs=None):
        lFiles = []
        if Dirs is None:
            for (dirpath, dirnames, filenames) in os.walk(self.home_dir[0]):
                lFiles.append(filenames)

        else:
            for (dirpath, dirnames, filenames) in os.walk(Dirs):
                lFiles.append(filenames)

        return lFiles



class jcFile:
    def getFileinfo(self, fpath=None):
		if fpath is not None:
			inName  = fpath.split('/');
			fPath   = inName[:(len(inName)-1)];
			fPath   = '/'.join(fPath);
			inName  = inName[(len(inName)-1)].split('.');
			fExt    = inName[(len(inName)-1)];
			inName  = inName[:(len(inName)-1)];
			fName   = '.'.join(inName);
			return fPath, fName, fExt



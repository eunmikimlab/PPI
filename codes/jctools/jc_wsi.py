"""
   Jong Cheol's Package contain various utilities 
"""
import os
import subprocess
import time
import re
import csv
import openslide
import json
import csv
import json
from PIL import ImageDraw
import jcutils


class jcImagePath:
    def __init__(self, fpath=None):
	self.fpath = fpath;
	
    '''
	This function reads coordinates downloaded from ImagePath DB 
    '''
    def readImagePathFile(self, fpath=None):
	self.fpath = fpath;
	slideID  = [];
	slideExt = [];
	ROI      = [];
	with open(fpath, 'rb') as f:
	    objCSV = csv.reader(f, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=False, escapechar='\\')
	    
	    '''
	    for tmp in objCSV:
		ROI.append(tmp[2]);
		slideID.append(tmp[10]);
		slideExt.append(tmp[18]);
	    
	    data = {'roi':ROI, 'id':slideID, 'format':slideExt}
	    '''
	    data = list(objCSV);
	return (data);

    '''
	This function draw ROIs with a specific resolution
    '''
    def drawROI(self, infile=None, outfile=None, coord=None, level=5, color='green', oformat='png'):
	#// read input file
	if infile is not None:
	    try:
		islide = openslide.open_slide(infile);
	    except OSError:
		print('cannot open', infile);
	else:
	    raise ValueError('Input file is not given');
	
	#// preparing for output file
	if outfile is None:
	    obj   = jcutils.jcFile();
	    p,n,e = obj.getFileinfo(infile);
	    outfile = '/'.join((p,n));
	    outfile = '.'.join((outfile,oformat));
	else:
	    obj     = jcutils.jcFile();
	    p,n,e   = obj.getFileinfo(outfile);
	    oformat = e;
	
	#// calculating downsample factors (http://openslide.org/api/python/)
	factor = islide.level_downsamples[level];
	
        #// build image with desired resolution
        img = islide.read_region((0,0),level,islide.level_dimensions[level])
	    
        #// drawing ROIs on an image
        imgID = ImageDraw.Draw(img)

	#// coordinate data must be dictionary data
	# - For polygon
	# coord = {'type':'polyline',
	#          'points': list...}
	#
	# - For circle
	# coord = {'type':'circle',
	#          'cx': <int>,
	#          'cy': <int>,
	#	   'r' : <int>}
	#
	# - For rect
	# coord = {'type':'rect',
	#          'x': <int>,
	#          'y': <int>,
	#	   'width' : <int>,
	#	   'height': <int>}
	
	if coord is not None:
	    #// for polygon
	    if coord['type'].lower() == 'polyline':
		PTs = coord['points'];
		for i in range(0,len(PTs)):
		    PTs[i] = int(PTs[i]/factor)

		imgID.polygon(PTs, outline=color);
		
	    #// for circle
	    if coord['type'].lower() == 'circle':
		cx = int(coord['cx']/factor);
		cy = int(coord['cy']/factor);
		r  = int(coord['r']/factor);
		imgID.ellipse((cx-r,cy-r,cx+r,cy+r), outline=color)
		
	    #// for rectangle
	    if coord['type'].lower() == 'rect':
		x = int(coord['x']/factor);
		y = int(coord['y']/factor);
		w = int(coord['width']/factor);
		h = int(coord['height']/factor);
		imgID.rectangle(((x, y), (x+w, y+h)), outline=color)
		
	img.save(outfile);
	islide.close();

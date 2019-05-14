#!/usr/bin/env python2

import os.path, subprocess, sys

def layer (fname, nsamples):
	ntimes = (nsamples/5000)-1
	fsplit = fname.split('-')
	fext = fsplit[3].split['.'][1]
	fout = fsplit[0] + fsplit[1] + str(nsamples) + fext
	with open(fname, "r") as fin:
		for line in fin:
			linecopy = str(line).rstrip('\n')
			fcopy = linecopy[linecopy.find('}')+1:len(linecopy)]
			for i in range(ntimes):
				linecopy += fcopy
			fout.write(linecopy + "\n")

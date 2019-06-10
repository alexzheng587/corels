#!/usr/bin/env python2
import os.path, os, subprocess, sys

from subprocess32 import STDOUT, check_output

import re

def natural_sort(l): 
	convert = lambda text: int(text) if text.isdigit() else text.lower() 
	alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key)  ] 
	return sorted(l, key = alphanum_key)

def run_corels (fname, thread_count, iterno):
	fargs = fname + ".out " + fname + ".label " + fname + ".minor"
	command = "../src/corels -c 2 -p 1 -r 0.01 -v 10 -n 2000000000 -i " + str(iterno) + " -t " + str(thread_count) + " " + fargs
	print(command)
	#o = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True).communicate()
	#return o[0].decode()
	output = check_output(command, stderr=STDOUT, timeout=21600, shell=True)
	return output.decode()

def main ():
	dirs = natural_sort(os.listdir('/data/corels/datasets/census/'))
	for i in dirs:
		if not i.endswith('.out'):
			continue
		fname = '/data/corels/datasets/census/' + i[:-4]
		print(fname)
		if not os.path.isfile(fname + '.out'):
			continue
		basename = fname[fname.rfind('/')+1:]
		for x in range(3):
			t = 1
			with open("../logs/proc-logs/" + basename + "-" + str(t) + "-" + str(x) + ".txt", "w") as f:
				f.write(run_corels(fname, t, x))

main()

#!/usr/bin/env python2
import os.path, subprocess, sys

def run_corels (fname, iterno):
	fargs = fname + ".out " + fname + ".label " + fname + ".minor"
	command = "~/bbverify/src/corels -c 2 -p 1 -r 0.01 -v progress,log -n 2000000000 -i " + str(iterno) + " " + fargs
	print(command)
	o = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True).communicate()
	return o[0].decode()

def main ():
	fname = "../data_check/census_c1s50ky0.3_10_train"
	basename = fname[fname.rfind('/')+1:]
	for x in range(1):
		with open("../logs/proc-logs/" + basename + "_" + str(x) + ".txt", "w") as f:
			f.write(run_corels(fname, x))

main()

#!/usr/bin/env python2
import os.path, os, subprocess, sys

def run_corels (fname, thread_count, iterno):
	fargs = fname + ".out " + fname + ".label " + fname + ".minor"
	command = "../src/corels -c 2 -p 1 -r 0.01 -v 10 -n 2000000000 -i " + str(iterno) + " -t " + str(thread_count) + " " + fargs
	print(command)
	o = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=None, shell=True).communicate()
	return o[0].decode()

def main ():
	for i in os.listdir('/home/asaligrama/census-data/20181224/use/'):
		if not i.endswith('.out'):
			continue
		fname = '/home/asaligrama/census-data/20181224/use/' + i[:-4]
		print(fname)
		if not os.path.isfile(fname + '.out'):
			continue
		basename = fname[fname.rfind('/')+1:]
		for t in [1, 2, 4, 8, 16, 32]:
			for x in range(5):
				with open("../logs/proc-logs/" + basename + "-" + str(t) + "-" + str(x) + ".txt", "w") as f:
					f.write(run_corels(fname, t, x))

main()

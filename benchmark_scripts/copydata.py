#!/usr/bin/env python
import sys

filename = sys.argv[1]

fsplit = filename.split('-')
card, rule, samp, ext = fsplit[0], fsplit[1], fsplit[2], fsplit[3]

with open(filename, 'r') as fin:
	for i in [2, 10, 20, 100, 200]:
		with open('-'.join([card, rule, str(i*int(samp)), ext]), 'w') as fout:
			fin.seek(0)
			for line in fin:
				loc = line.index('}')+1
				lout = line[:loc]
				for j in range(i):
					lout += line[loc:-1]
				fout.write(lout + '\n')
				fout.flush()
			fout.close()
	fin.close()

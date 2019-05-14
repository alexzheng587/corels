#!/usr/bin/env python3

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
import numpy as np
import csv

fig = plt.figure()
ax = plt.axes(projection='3d')

X = []
Y = []
Z = []

with open('..\data\scalability.csv') as fin:
	reader = csv.reader(fin)
	next(reader, None)
	for line in reader:
		xs = float(line[1])
		ys = float(line[0])
		zs = float(line[2])
		print(xs, ys, zs)
		X.append(xs)
		Y.append(ys)
		Z.append(zs)

X, Y, Z = np.array(X), np.array(Y), np.array(Z)

surf = ax.plot_trisurf(np.log10(X), np.log10(Y), np.log10(Z), cmap=cm.coolwarm, linewidth=0, antialiased=True)

ax.set_xlabel('Dataset size log10(#rules)')
ax.set_ylabel('Dataset size log10(#samples)')
ax.set_zlabel('Runtime log10(s)')

fig.colorbar(surf, shrink=0.5, aspect=5)

plt.show()

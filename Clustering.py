# needed imports
import sys
import os
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
import numpy as np
from scipy.cluster.hierarchy import fcluster

#IMPORTAR COORDS en X

Coord=[]
Dens=[]

STRUCTUREPATH = os.path.abspath(sys.argv[1])
STRUCTURENAME = STRUCTUREPATH.split('/')[-1].split('.')[0]
filein=open(STRUCTUREPATH, "r")
lines=filein.readlines()
filein.close()

for i in lines:
	if "CP" in i:
		k=lines.index(i)
		k+=1
		while lines[k]!="\n":
			split=lines[k].split()
			Coord.append([float(split[1]), float(split[2]), float(split[3])])
			Dens.append(float(split[4]))
			k+=1
#

Y=np.array(Coord)


Z=linkage(Y, 'ward')
max_d=1
clusters = fcluster(Z, max_d, criterion='distance')

#print clusters
#print max(clusters)

CP=[]

ncp=[[]]*max(clusters)



j=0
for i in clusters:
	n=i-1
	ncp[n]=ncp[n]+[[Coord[j][0],Coord[j][1],Coord[j][2],Dens[j]]]
	j=j+1


ccp=[]


for i in ncp:
	x=[]
	y=[]
	z=[]
	d=[]
	for j in i:
		x.append(j[0])
		y.append(j[1])
		z.append(j[2])
		d.append(j[3])
	X=sum(x)/len(x)
	Y=sum(y)/len(y)
	Z=sum(z)/len(z)
	D=sum(d)
	ccp.append([X, Y, Z, D])

with open(os.path.abspath(os.getcwd())+"/SPHERES_CP.py", "w") as out:
	out.write("from pymol.cgo import *\nfrom pymol import cmd\n\nspherelist = [\n")
	for i in ccp:
		out.write("\tCOLOR,\t0.100,\t1.000,\t0.000,\n\tSPHERE, \t%6.3f,\t%6.3f,\t%6.3f,\t%6.3f,\n"%(i[0], i[1], i[2], i[3]))
	out.write("\t]\n\ncmd.load_cgo(spherelist, 'segment',   1)")

out.close()

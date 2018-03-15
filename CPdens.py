import sys, getopt
import os
import matplotlib.pyplot as plt
plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
from subprocess import Popen, PIPE

from metilar import metilar
from Com_gen import getCOM
from Densidadtotal import densCP
from symbol_to_vdw import symbol_to_vdw

#############################################
# Reads the inputfile and the stage selected
#############################################
argv= sys.argv[1:]

inputfile = ''
stage = ''
try:
	opts, args = getopt.getopt(argv,"hi:s:",["ifile=","s="])
except getopt.GetoptError:
	print 'CPdens.py -i <inputfile> -s <stage>'
	sys.exit(2)
for opt, arg in opts:
	if opt == '-h':
		print 'CPdens.py -i <inputfile> -s <stage>'
		sys.exit()
	elif opt in ("-s", "--ifile"):
		stage = arg
	elif opt in ("-i", "--stage"):
		inputfile = arg

filein = open(os.path.abspath(os.getcwd())+"/"+inputfile,"r")
lines = filein.readlines()
filein.close()

structure=lines[0].split()[1]
ligand=lines[1].split()[1]
residues=[]
[residues.append(float(i)) for i in lines[2].split()[1:]]
optmethod=lines[3].split()[1]
g09method=lines[4].split()[1]
charge=lines[5].split()[1]
multS=lines[6].split()[1]
nprocs=lines[7].split()[1]

resnames=[]
resdens=[]
frames=[]
dtot=[]

#############################################
# Generate .mop and .aux files from pdb structures
#############################################
if stage=="1":
	
	[frames.append(i.split(".")[0]) for i in os.listdir(os.getcwd()) if (structure in i and ".pdb" in i)]

	for i in frames:
		metilar(i, optmethod, residues, ligand, charge)
		
#############################################
# Generate .com files from mopac output
#############################################
elif stage=="2":
	
	[frames.append(i.split(".")[0]) for i in os.listdir(os.getcwd()) if (structure in i and ".out" in i)]
	
	for i in frames:
		getCOM(i, g09method, charge, multS, nprocs)
		
#############################################
# Generate .txt files from .wfn
#############################################
elif stage=="3":
	
	[frames.append(i.split(".")[0]) for i in os.listdir(os.getcwd()) if (structure in i and ".wfn" in i)]
	for i in frames:
		p = Popen(["/home/esteban/Sulfonamidas_sulfamatos/AIM/hCAII/programa/Multiwfn", i+".wfn"], stdout=PIPE,stdin=PIPE, stderr=PIPE)
		print p.communicate("2\n3\n7\n-1")
		os.rename("CPprop.txt", i+".txt")
		
#############################################
# Generate a unique CPs log file from all frames
#############################################
elif stage=="4":
	densities=[]
	[frames.append(i.split(".")[0]) for i in os.listdir(os.getcwd()) if (structure in i and ".txt" in i)]
	with open(os.path.abspath(os.getcwd())+"/"+structure+"_CP.log", "w") as out:
		out.write("Ligand-Protein bond critical point analysis\n"+"-"*40+"\n")
	out.close()
	
	n=0
	for i in frames:
		m=0
		dat=densCP(i)
		if n==0:
			resnames=dat[1]
			resdens.append(dat[2])
			dtot.append(dat[0])
		else:
			for j in dat[2]:
				resdens.append(dat[2])
				m+=1
			dtot.append(dat[0])
		n+=1
	
	
	y_pos=np.array(len(resnames))
	mean_resdens=[]
	std_resdens=[]
	
	for i in range(len(resnames)):
		a=[]
		for j in range(len(dtot)):
			a.append(resdens[j][i])
		mean_resdens.append(np.mean(a))
		std_resdens.append(np.std(a))
	

	density=np.mean(dtot)
	density_std=np.std(dtot)

	with open(os.path.abspath(os.getcwd())+"/"+structure+"_CP.log", "a") as out:
		out.write(("-"*40+"\nMean total density: %7.4f +- %7.4f\n"+"-"*40+"\n")%(density, density_std))
	
	print("Se creo el archivo "+structure+"_CP.log")
		
	#Grafica la decomposicion

	plt.rcdefaults()
	fig, ax = plt.subplots()
    
	 
	y_pos = np.arange(len(resnames)-1)
	y_val = np.array(mean_resdens[0:-1])
	err=np.array(std_resdens[0:-1])
    
	ax.barh(y_pos, y_val, align='center',
			color='green', ecolor='black', xerr=err)
	ax.set_yticks(y_pos)
	ax.set_yticklabels(resnames)
	ax.invert_yaxis()  # labels read top-to-bottom
	ax.set_xlabel('Densidad')
	ax.set_title('Decomposicion por residuo de la densidad electronica\nDensidad total = '+str(density))
	ax.set_xlim(left=0, right=0.12)
    
	pylab.savefig(structure+'.png')
		
	print("Se creo el archivo "+structure+".png")
	


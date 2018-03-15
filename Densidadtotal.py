#### Funcion para extraer-filtrar-clasificar puntos criticos ####
# 1. Utiliza el archivo auxiliar para asociar los atomos con sus 
#	 respectivos residuos (instanciados dentro de la clase residue)
# 2. Almacena las coordenadas y la densidades de todos los puntos 
#    criticos del tipo (3 -1)
# 3. Encuentra la transformacion de cooords del Gaussian
# 4. Asocia los PC con los residuos y suma las densidades
# 5. Escribe en la salida
# 6. Devuelve el valor de Dtot y las listas de residuos/densidades

#Ultima modificacion: 16-02-18

import sys
import os
import matplotlib.pyplot as plt
plt.rcdefaults()
import numpy as np
import matplotlib.pyplot as plt
import pylab
import math
from symbol_to_vdw import symbol_to_vdw

def densCP(frame):


	NAME=frame.split("_")[0]+"_"+frame.split("_")[1]
	


	ATOMS = {"1":"H ", #La salida del Gaussian usa numeros en vez de simbolo
	"6":"C ",
	"7":"N ",
	"8":"O ",
	"30":"Zn",
	"16":"S "}
	LOG = []
	R = []
	reslist= []
	Dtot= 0
	CPS_a = []
	CPS = []
	CENT = []

	#Define la clase donde se contendran los atomos y la densidad electonica

	class residue:
		def __init__(self, resname, resid):
			self.resname = resname
			self.resid = resid
			self.Dpart = 0
			self.Atomlst = []
		def contain(self, atmline):
			if atmline.split()[2]==self.resid:
				self.Atomlst.append(atmline.split()[0])
		

	#Instancia los residuos a partir del archivo auxiliar

	auxfile = open(os.path.abspath(os.getcwd())+"/"+frame+".aux","r")
	lines = auxfile.readlines()
	auxfile.close()

	for i in lines:
		 if i.split()[2] not in reslist:
			 reslist.append(i.split()[2])
			 R.append(residue(i.split()[1],i.split()[2]))


	for i in lines:
		for j in R:
			j.contain(i)
			

	#-----------------------------------------------------------------------
	#Abre el CP, suma las densidades y guarda las coord

	cpfile = open(os.path.abspath(os.getcwd())+"/"+frame+".txt","r")
	lines = cpfile.readlines()
	cpfile.close()


	for i in lines:
		if "================   CP" in i and "(3,-1)" in i:
			indexcoord=lines.index(i)
			indexcoord +=2
			dens=float(lines[indexcoord].split()[4])
			indexcoord -=1
			xcoord=float(lines[indexcoord].split()[2])*0.529
			ycoord=float(lines[indexcoord].split()[3])*0.529
			zcoord=float(lines[indexcoord].split()[4])*0.529
			
			
			CPS_a.append({"x":xcoord, "y":ycoord, "z":zcoord,
			"dens":dens})
			

		
	#-----------------------------------------------------------------------
	#Abre el .log de gaussian y encuentra la transformacion de coords
	
	cpfile = open(os.path.abspath(os.getcwd())+"/"+frame+".log","r")
	lines = cpfile.readlines()
	cpfile.close()
	
	indexcoord=lines.index(" Symbolic Z-matrix:\n")
	indexcoord += 2
	i0=indexcoord
	i1=indexcoord+1
	i2=indexcoord+2
	i3=indexcoord+3
	
	r1 = np.array([float(lines[i0].split()[1]), float(lines[i0].split()[2]), float(lines[i0].split()[3])])
	r2 = np.array([float(lines[i1].split()[1]), float(lines[i1].split()[2]), float(lines[i1].split()[3])])
	r3 = np.array([float(lines[i2].split()[1]), float(lines[i2].split()[2]), float(lines[i2].split()[3])])
	r4 = np.array([float(lines[i3].split()[1]), float(lines[i3].split()[2]), float(lines[i3].split()[3])])
	
	indexcoord=lines.index(" Number     Number       Type             X           Y           Z\n")
	indexcoord += 2
	i0=indexcoord
	i1=indexcoord+1
	i2=indexcoord+2
	i3=indexcoord+3
	
	r1t = np.array([float(lines[i0].split()[3]), float(lines[i0].split()[4]), float(lines[i0].split()[5])])
	r2t = np.array([float(lines[i1].split()[3]), float(lines[i1].split()[4]), float(lines[i1].split()[5])])
	r3t = np.array([float(lines[i2].split()[3]), float(lines[i2].split()[4]), float(lines[i2].split()[5])])
	r4t = np.array([float(lines[i3].split()[3]), float(lines[i3].split()[4]), float(lines[i3].split()[5])])
	
	r1_0 = np.array([0, 0, 0])
	r2_0 = r2 - r1
	r3_0 = r3 - r1
	r4_0 = r4 - r1
	
	r1t_0 =np.array([0,0,0])
	r2t_0 = r2t - r1t
	r3t_0 = r3t - r1t
	r4t_0 = r4t - r1t
	
	Q = np.stack((r2_0, r3_0, r4_0), axis=0)
	P = np.stack((r2t_0, r3t_0, r4t_0), axis=0)
	M = np.linalg.solve(P, Q)
	
	
	#Guarda las coordenadas del log de gaussian

	while lines[indexcoord]!=" ---------------------------------------------------------------------\n":
		sym=ATOMS[lines[indexcoord].split()[1]]
		xcoord=float(lines[indexcoord].split()[3])
		ycoord=float(lines[indexcoord].split()[4])
		zcoord=float(lines[indexcoord].split()[5])
		LOG.append({"sym":sym, "x":xcoord, "y":ycoord, "z":zcoord})
		indexcoord +=1

	#Busca los pares a R<suma_rVdW
	#Guarda las coord de los puntos medios entre atomos de la prot y del lig

	primer=int(R[-1].Atomlst[0])
	ultimo=int(R[-1].Atomlst[-1])
	a=0
	for i in LOG[0:(primer-1)]: #Itera sobre las coord de todos los atomos menos el ultimo residuo
		a +=1
		b=0
		for j in LOG[(primer-1):ultimo]: #Itera sobre las coord de los atomos del ultimo residuo (ligando)
			b+=1
			r= math.sqrt((i["x"]-j["x"])**2+(i["y"]-j["y"])**2+(i["z"]-j["z"])**2)
			vdw=symbol_to_vdw(i["sym"].strip())+symbol_to_vdw(j["sym"].strip())
			
			if r<=(vdw*1.15):	#Si el factor de tolerancia es grande pueden entrar puntos criticos extras
				factor=symbol_to_vdw(j["sym"].strip())/vdw
				xcent=(i["x"]-j["x"])*factor+j["x"]
				ycent=(i["y"]-j["y"])*factor+j["y"]
				zcent=(i["z"]-j["z"])*factor+j["z"]
				CENT.append({"x":xcent, "y":ycent, "z":zcent, "idx":a, "idl":b}) 

	#Asocia los puntos criticos con el centro mas cercano y filtra una distancia maxima

	for i in CPS_a:
		dmin=10
		a=0
		for j in CENT:
			d=math.sqrt((i["x"]-j["x"])**2+(i["y"]-j["y"])**2+(i["z"]-j["z"])**2)
			if d < dmin:
				dmin=d
				i["idx"]=j["idx"]
				i["idl"]=j["idl"]

		if dmin < 0.4: #Distancia maxima permitida entre centro y PC
			CPS.append(i)
			
	#Suma las densidades al residuo que le corresponda
	
	for i in CPS:
		for j in R:
			if str(i["idx"]) in j.Atomlst:
				j.Dpart+=i["dens"]
				
	#Suma la densidad total
	
	Dtot=0
	for j in R:
		Dtot+=j.Dpart
		
	for i in CPS:
		F=np.array([i["x"], i["y"], i["z"]])
		I=np.dot((F-r1t), M) + r1
		i["x"]=I[0]
		i["y"]=I[1]
		i["z"]=I[2]
		
		
	#Escribe en el archivo de salida la informacion del frame
	
	with open(os.path.abspath(os.getcwd())+"/"+NAME+"_CP.log", "a") as out:
		n=1
		out.write(("-"*40+"\nFrame\t%s\n"+"-"*40+"\nTotal Density: %7.4f\n")%(frame.split("_")[2], Dtot))
		out.write("Residue\t\tDensity\n")
		for i in R:
			out.write("%s %s\t\t%7.4f\n"%(i.resname, i.resid, i.Dpart))
		out.write("CP\t\tX\t\tY\t\tZ\t\tDens\n")
		for i in CPS:
			out.write("%i\t%7.4f\t%7.4f\t%7.4f\t%7.4f\t%i\t%i\n"%(n, i["x"], i["y"], i["z"], i["dens"], i["idx"], i["idl"]))
			n+=1
		out.write("\n")
	out.close()
	
	#Devuelve el valor de densidad total y listas de residuos/densidades
	
	resnames = []
	densities = []
	[resnames.append(i.resname+i.resid) for i in R]
	[densities.append(i.Dpart) for i in R]
	return [Dtot, resnames, densities]
    

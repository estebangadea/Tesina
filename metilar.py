#Funcion para cambiar los atomos del backbone de los residuos por metilos

import sys
import os

def metilar(frame, method, residues, ligand, charge):

	PDBPATH = os.path.abspath(os.getcwd())+"/"+frame+".pdb"
	atoms=[]
	CAcoord=[0]*260
	
	#Abre el pdb y guarda las lineas
	
	pdbfile = open(PDBPATH,"r")
	lines = pdbfile.readlines()
	pdbfile.close()
	
	i=0
	while i <= len(lines): 
		if ("ATOM" not in lines[i] and "HETATM" not in lines[i]):
			del lines[i]
		i+=1
	
	#Itero sobre las lineas guardando elemento y coordenadas
	for i in lines:
		if int(i.split()[5]) in residues and i.split()[2]=="CA":
			xca=float(i.split()[6])
			yca=float(i.split()[7])
			zca=float(i.split()[8])
			res=int(i.split()[5])
			CAcoord[res]={"x":xca, "y":yca, "z":zca}
	for i in lines:
		if int(i.split()[5]) in residues:
			xcoord=float(i.split()[6])
			ycoord=float(i.split()[7])
			zcoord=float(i.split()[8])
			at=i.split()[2]
			at=at[0]
			if at=="Z":
				at="ZN"
			res=int(i.split()[5])
			resn=i.split()[3]
			if res-1 in residues and res+1 in residues: #Si el residuo tiene su anterior y su posterior hay que transcribirlo completo
				atoms.append({"x":xcoord, "y":ycoord, "z":zcoord, "at":at, "res":res, "resn":resn})
				
			elif res-1 in residues: #Si esta el residuo anterior solo quiero quitar el O y reemplazar el carbono 
				if i.split()[2]=="C" and i.split()[3]!=ligand:
					at="H"
					cax=CAcoord[res]["x"]
					cay=CAcoord[res]["y"]
					caz=CAcoord[res]["z"]
					xcoord= (xcoord-cax)*1.09/1.54+cax
					ycoord= (ycoord-cay)*1.09/1.54+cay
					zcoord= (zcoord-caz)*1.09/1.54+caz
					atoms.append({"x":xcoord, "y":ycoord, "z":zcoord, "at":at, "res":res, "resn":resn})
				elif i.split()[2]=="O" and i.split()[3]!=ligand:
					j=None
				else:
					atoms.append({"x":xcoord, "y":ycoord, "z":zcoord, "at":at, "res":res, "resn":resn})
					
			elif res+1 in residues: #Si esta el residuo posterior quiero quitar el H y reemplazar el H
				if i.split()[2]=="H" and i.split()[3]!=ligand:
					j=None
				elif i.split()[2]=="N" and i.split()[3]!=ligand:
					at="H"
					cax=CAcoord[res]["x"]
					cay=CAcoord[res]["y"]
					caz=CAcoord[res]["z"]
					xcoord= (xcoord-cax)*1.09/1.47+cax
					ycoord= (ycoord-cay)*1.09/1.47+cay
					zcoord= (zcoord-caz)*1.09/1.47+caz
					atoms.append({"x":xcoord, "y":ycoord, "z":zcoord, "at":at, "res":res, "resn":resn})
				else:
					atoms.append({"x":xcoord, "y":ycoord, "z":zcoord, "at":at, "res":res, "resn":resn})
			else:  #En caso que no haya residuos contiguos hay que metilarlo completo
				if i.split()[2]=="C" and i.split()[3]!=ligand:
					at="H"
					cax=CAcoord[res]["x"]
					cay=CAcoord[res]["y"]
					caz=CAcoord[res]["z"]
					xcoord= (xcoord-cax)*1.09/1.54+cax
					ycoord= (ycoord-cay)*1.09/1.54+cay
					zcoord= (zcoord-caz)*1.09/1.54+caz
					atoms.append({"x":xcoord, "y":ycoord, "z":zcoord, "at":at, "res":res, "resn":resn})
				elif i.split()[2]=="O" and i.split()[3]!=ligand:
					j=None
				elif i.split()[2]=="H" and i.split()[3]!=ligand:
					j=None
				elif i.split()[2]=="N" and i.split()[3]!=ligand:
					at="H"
					cax=CAcoord[res]["x"]
					cay=CAcoord[res]["y"]
					caz=CAcoord[res]["z"]
					xcoord= (xcoord-cax)*1.09/1.47+cax
					ycoord= (ycoord-cay)*1.09/1.47+cay
					zcoord= (zcoord-caz)*1.09/1.47+caz
					atoms.append({"x":xcoord, "y":ycoord, "z":zcoord, "at":at, "res":res, "resn":resn})
				else:
					atoms.append({"x":xcoord, "y":ycoord, "z":zcoord, "at":at, "res":res, "resn":resn})
	
	
	#Escribo un archivo auxiliar para la separacion por residuos
	n=1
	with open(os.path.abspath(os.getcwd())+"/"+frame+".aux", "w") as out:
		for i in atoms:
			out.write("%i \t %s \t %i \t %s\n"%(n, i["resn"], i["res"], i["at"]))
			n+=1
	out.close()	
	
	#Escribo la salida para el mopac
	with open(os.path.abspath(os.getcwd())+"/"+frame+".mop", "w") as out:
		out.write(method+" NOOPT OPT-H CHARGE="+charge+"\n\n\n")
		for i in atoms:
			out.write("%s \t%7.4f \t%7.4f \t%7.4f \n"%(i["at"], i["x"], i["y"], i["z"]))
	out.close()
	print "Se creo el archivo "+frame+".mop"

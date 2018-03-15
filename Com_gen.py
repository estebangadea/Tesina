#FUncion para concatenar trbajos de MOPAC con Gaussian
#Extrae las coordenadas del .out y escribe el nuevo com


import sys
import os

def getCOM(frame, method, charge, multS, nprocs):

	OUTPATH = os.path.abspath(os.getcwd())+"/"+frame+".out"

	#1. Abre el .out y almacena las coordenadas cartesianas finales
	atomlst=[]
	try:
		mopacfile = open(OUTPATH,"r")
		lines = mopacfile.readlines()
		mopacfile.close()
		try:
			indexcoord=lines.index("                             CARTESIAN COORDINATES\n")
			indexcoord += 2
			while lines[indexcoord]!="\n":
				linediv = lines[indexcoord].split()
				symbol = linediv[1]
				number = int(linediv[0])
				xcor = float(linediv[2])
				ycor = float(linediv[3])
				zcor = float(linediv[4])
				atomlst.append({"sym":symbol, "num":number, "x":xcor, "y":ycor, "z":zcor})
				indexcoord += 1
			#2. Escribe el nuevo .com
			with open(os.path.abspath(os.getcwd())+"/"+frame+".com", "w") as out:
				out.write("%mem=100MW\n%nprocs="+nprocs+"\n#n "+method+" SP OUTPUT=WFN SCF=MAXCYCLE=1000\n\nC-alfas a metilos opt-H\n\n"+multS+" "+charge+"\n")
				for i in atomlst:
					out.write("%s \t%10.5f \t%10.5f \t %10.5f\n"%(i["sym"],i["x"],i["y"],i["z"]))
				out.write("\n%s.wfn\n\n"%(frame))
			print("Se creo el archivo "+frame+".com")
		except(ValueError):
			print ("Parece que el archivo "+frame+".out no se pudo leer correctamente")
	except (IOError, OSError):
		print ("No se pudo abrir el archivo "+frame+".out")


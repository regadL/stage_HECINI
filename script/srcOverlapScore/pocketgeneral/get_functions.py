#import os
#import string
#import sys
#import re
#from math import *
#from PDB import *
#from random import *
#from numpy import *
#from Scientific.Functions.LeastSquares import leastSquaresFit
#from Scientific.Geometry.Objects3D import *
#from subprocess  import Popen
#
#import pathDirectory
#import runOtherProg
#import parseNACCESS
#
#extend = 5
#
#
#
#
#
#
#
#
#
#def angle(vpa,vpb):
#	a0 = vpa[0]
#	a1 = vpa[1]
#	a2 = vpa[2]
#	b0 = vpb[0]
#	b1 = vpb[1]
#	b2 = vpb[2]
#	scalaire = a0*b0+a1*b1+a2*b2
#	# voir dot ... ?
#	return arccos(scalaire)
#
#def mini2(a,b):
#	if a<=b:
#		return a
#	return b
#
#def mini3(a,b,c):
#	if a<=b:
#		return mini2(a,c)
#	return mini2(b,c)
#
#
#def maxi2(a,b):
#	if a<=b:
#		return b
#	return a
#
#def maxi3(a,b,c):
#	if a<=b:
#		return maxi2(b,c)
#	return maxi2(a,c)
#	
#
#
#
## calcul des composantes principales pour le protomol et pour le ligand
## separe dans get_descriptors en fait
##def principal_components(fichier,pock,fileout,proto):
##	#compute_lambdas(fichier,pock,'pocket',"protomol/"+proto+"_"+pock+fichier+"-residues_around_4_0.pdb",fileout)
##	compute_lambdas(fichier,pock,'pocket',"protomol/"+proto+"_"+pock+"-protomol.pdb",fileout)
##	if (os.path.isfile("ligand.pdb")!=0):
##		compute_lambdas(fichier,pock,'ligand',"ligand.pdb",fileout)
##
##	
#	
##def donors(seq):
##	"""
##	Retrieve number of h bond donor, base sequences
##	args: -> seq
##	return: -> integer count
##	"""
##	nb_d = 0
##	for aa_d in HBd:
##		nb_d += string.count(seq,aa_d)
##	print nb_d, seq
##	return nb_d
#
#
#
## calcul du polarity ratio - mis de cote ... a mon avis il manque les C.2 !!
## et puis on veut faire sur le ligand et la poche en terme de residus donc il faut tenir compte des S
## definition du polarity ratio : Transient Pocket ... Eyrisch and Helms - J.Med.Chem. 2007
##def polarity_ratio_old(flag,fileout):
##	os.system("grep \"O.2\" protomol/protomol.mol2 | wc -l > temp_O.out")
##	fileO = open("temp_O.out",'r')
##	nO = fileO.readline()
##	nO = int(string.replace(nO,'\n',''))
##	fileO.close()
##	os.system("rm temp_O.out")
##	os.system("grep \"N.3\" protomol/protomol.mol2 | wc -l > temp_N.out")
##	fileN = open("temp_N.out",'r')
##	nN = fileN.readline()
##	nN = int(string.replace(nN,'\n',''))
##	fileN.close()
##	os.system("rm temp_N.out")
##	os.system("grep \"C.3\" protomol/protomol.mol2 | wc -l > temp_C.out")
##	fileC = open("temp_C.out",'r')
##	nC = fileC.readline()
##	nC = int(string.replace(nC,'\n',''))
##	fileC.close()
##	os.system("rm temp_C.out")
##	ratio = (nO+nN)/float(nO+nN+nC)
##	ratio = "%.2f" %ratio
##	fileout.write(flag+"_polarity_ratio_old\t"+ratio+"\n")
#
## manque les C.2 dans le precedant ... nouveau calcul .. en fait C.2 = O donc je sais pas
##def ratio(path_file_pocket, path_result, type_study, filout):
##	
##	cmd_O = "grep \" O \" "+path_file_pocket+" | wc -l > " + path_result + "temp_O.out"
##	os.system(cmd_O)
##	fileO = open(path_result + "temp_O.out",'r')
##	nO = fileO.readline()
##	nO = int(string.replace(nO,'\n',''))
##	fileO.close()
##	os.system("rm " + path_result + "temp_O.out")
##	
##	os.system("grep \" N \" "+path_file_pocket+" | wc -l > "  + path_result + "temp_N.out")
##	fileN = open( path_result +"temp_N.out",'r')
##	nN = fileN.readline()
##	nN = int(string.replace(nN,'\n',''))
##	fileN.close()
##	os.system("rm "+ path_result + "temp_N.out")
##	
##	os.system("grep \" S \" " + path_file_pocket+ " | wc -l > "+ path_result + "temp_S_all.out")
##	# sans les S des methionine ... cf Olivier
##	os.system("grep \" S..*CYS\" "+path_file_pocket+" | wc -l > "+ path_result + "temp_S.out")
##	fileS = open( path_result + "temp_S.out",'r')
##	nS = fileS.readline()
##	nS = int(string.replace(nS,'\n',''))
##	fileS.close()
##	os.system("rm " + path_result + "temp_S.out")
##	fileSall = open(path_result +"temp_S_all.out",'r')
##	nSall = fileSall.readline()
##	nSall = int(string.replace(nSall,'\n',''))
##	fileSall.close()
##	os.system("rm " + path_result + "temp_S_all.out")
##	os.system("grep \" C..\" "+path_file_pocket+" | wc -l > "+ path_result + "temp_C.out")
##	fileC = open(path_result + "temp_C.out",'r')
##	nC = fileC.readline()
##	nC = int(string.replace(nC,'\n',''))
##	fileC.close()
##	os.system("rm " +path_result + "temp_C.out")
##	
##	# polarity ratio
##	print nO, nN, nSall, nC
##	
##	ratio = (nO+nN+nS)/float(nO+nN+nSall+nC)
##	# hydrophobicity ratio voir Burgoyne et Jackson Bioinformatics
##	ratio2 = (nC+nS)/float(nO+nN+nSall+nC)
##	#ratio = "%.2f" %ratio
##	ratio = str(ratio)
##	ratio2 = "%.2f" %ratio2
##	
##	filout.write("pocket_polarity_ratio_" + type_study + "\t"+ratio+"\n")
##	filout.write("pocket_hydrophobicity_ratio_"+ type_study + "\t"+ratio2+"\n")
#
#
#def ratioOriginal(fichierin, flag, filout, na = 0):
#	if na==0:
#		point = 1
#		if flag == "protomol":
#			fichier = fichierin
#		elif flag == "ligand":
#			fichier = "ligand"+fichierin+".mol"
#			point = 0
#		elif flag == "pocket":
#			fichier = runOtherProg.babelPDBtoMOL2(fichierin)
#			
#		elif flag == "protein":
#			fichier = "protein.mol2"
#		else:
#			print "flag inconnu"
#		if point == 1:
#			os.system("grep \" O\.\" "+fichier+" | wc -l > temp_O.out")
#			fileO = open("temp_O.out",'r')
#			nO = fileO.readline()
#			nO = int(string.replace(nO,'\n',''))
#			fileO.close()
#			os.system("rm temp_O.out")
#			os.system("grep \" N\.\" "+fichier+" | wc -l > temp_N.out")
#			fileN = open("temp_N.out",'r')
#			nN = fileN.readline()
#			nN = int(string.replace(nN,'\n',''))
#			fileN.close()
#			os.system("rm temp_N.out")
#			os.system("grep \" S\.\" "+fichier+" | wc -l > temp_S_all.out")
#			# sans les S des methionine ... cf Olivier
#			os.system("grep \" S\..*CYS\" "+fichier+" | wc -l > temp_S.out")
#			fileS = open("temp_S.out",'r')
#			nS = fileS.readline()
#			nS = int(string.replace(nS,'\n',''))
#			fileS.close()
#			os.system("rm temp_S.out")
#			fileSall = open("temp_S_all.out",'r')
#			nSall = fileSall.readline()
#			nSall = int(string.replace(nSall,'\n',''))
#			fileSall.close()
#			os.system("rm temp_S_all.out")
#			os.system("grep \" C\.3\" "+fichier+" | wc -l > temp_C.out")
#			fileC = open("temp_C.out",'r')
#			nC = fileC.readline()
#			nC = int(string.replace(nC,'\n',''))
#			fileC.close()
#			os.system("rm temp_C.out")
#		elif point == 0:
#			os.system("grep \" O \" "+fichier+" | wc -l > temp_O.out")
#			fileO = open("temp_O.out",'r')
#			nO = fileO.readline()
#			nO = int(string.replace(nO,'\n',''))
#			fileO.close()
#			os.system("rm temp_O.out")
#			os.system("grep \" N \" "+fichier+" | wc -l > temp_N.out")
#			fileN = open("temp_N.out",'r')
#			nN = fileN.readline()
#			nN = int(string.replace(nN,'\n',''))
#			fileN.close()
#			os.system("rm temp_N.out")
#			os.system("grep \" S \" "+fichier+" | wc -l > temp_S.out")
#			fileS = open("temp_S.out",'r')
#			nS = fileS.readline()
#			nS = int(string.replace(nS,'\n',''))
#			nSall = nS
#			fileS.close()
#			os.system("rm temp_S.out")
#			os.system("grep \" C \" "+fichier+" | wc -l > temp_C.out")
#			fileC = open("temp_C.out",'r')
#			nC = fileC.readline()
#			nC = int(string.replace(nC,'\n',''))
#			fileC.close()
#			os.system("rm temp_C.out")
#		# polarity ratio -> a revoir plus propre
#		try : 
#			ratio = (nO+nN+nS)/float(nO+nN+nSall+nC)
#			ratio = "%.2f" %ratio
#		except :
#			ratio = "NA"
#		# hydrophobicity ratio voir Burgoyne et Jackson Bioinformatics
#		try :
#			ratio2 = (nC+nS)/float(nO+nN+nSall+nC)
#			ratio2 = "%.2f" %ratio2
#		except : 
#			ratio2 = "NA"
#	elif na==1:
#		ratio = "NA"
#	filout.write("pocket_polarity_ratio_" + flag + "\t"+ratio+"\n")
#	filout.write("pocket_hydrophobicity_ratio_" + flag + "\t"+ratio2+"\n")
#
## autre methode de recuperation des atomes accessibles ... je pense mieux mais de cote pour le moment car je n'ai pas envie de tout refaire
## maintenant
## Si on choisit cette methode un jour ... recoder pour le ligand
##def ratio_new(filout):
##	without = PDB("naccess/protein.asa")
##	within = PDB("protomol/protein_protomol_merged.asa")
##	acc_out = []
##	atom = []
##	acc_in = []
##	# accessibilite sans les sondes protomols
##	for resout in without:
##		for atout in resout:
##			if atout.header() == "ATOM":
##				temp = string.split(atout.occ())
##				acc_out.append(float(temp[0]))
##				atom.append(atout.atmName()[0])
##	# accessibilite avec les sondes protomol
##	for resin in within:
##		for atin in resin:
##			if atin.header() == "ATOM":
##				temp = string.split(atin.occ())
##				acc_in.append(float(temp[0]))
##	# difference
##	acc_diff = acc_out
##	for i in range(len(acc_out)):
##		acc_diff[i] = acc_out[i]-acc_in[i]
##	# compteur des atomes dont l'accessibilite varie
##	n = 0
##	# compteur des carbon et sulphur
##	ncs = 0
##	nnos = 0
##	for i in range(len(acc_diff)):
##		if acc_diff[i] != 0.:
##			if atom[i] != 'H':
##				n += 1
##			if atom[i] == 'C' or atom[i] == 'S':
##				ncs += 1
##			if atom[i] == 'N' or atom[i] == 'O' or atom[i] == 'S':
##				nnos += 1
##	ratio = ncs/float(n)
##	ratio2 = nnos/float(n)
##	ratio = "%.2f" %ratio
##	ratio2 = "%.2f" %ratio2
##	filout.write("pocket_hydrophobicity_ratio\t"+ratio+"\n")
##	filout.write("pocket_polarity_ratio\t"+ratio2+"\n")
#
#
#def dist(xa,ya,za,xb,yb,zb):
#	return sqrt((xa-xb)**2+(ya-yb)**2+(za-zb)**2)
#
## FONCTIONS pour la rugosite
## recuperation des numeros des atomes a 5 angs de l'atome numerote atomei
#def get_atoms_around(atomei,path_file_pdb):
#	ex = 5.
#	liste_atomes = []
#	prot = PDB(path_file_pdb)
#	for res in prot:
#		for at in res:
#			if int(at.atmNum()) == int(atomei):
#				xyz = at.xyz()
#				xi = xyz[0]
#				yi = xyz[1]
#				zi = xyz[2]
#				break
#	for res in prot:
#		for at in res:
#			xyzc = at.xyz()
#			xc = xyzc[0]
#			yc = xyzc[1]
#			zc = xyzc[2]
#			dist = sqrt((xc-xi)**2+(yc-yi)**2+(zc-zi)**2)
#			a = int(at.atmNum())
#			if at.header() == 'ATOM':
#				if (a != int(atomei) and dist <= ex):
#					liste_atomes.append(a)
#	return liste_atomes
#
#
#def getAtomPocketACCForArea(pocket_acc_parsed, path_file_asa, path_result, debug = 1):
#	
#	asa_parsed = parseNACCESS.fileASA(path_file_asa)
#	
#	liste_residus_pocket = []
#	list_atom_retrieve = []
#	
#	# on recupere les residus accessibles
#	for res in pocket_acc_parsed:
#		num_res = int(res.rNum())
#		liste_residus_pocket.append(num_res)
#	if debug : 
#		print liste_residus_pocket, "list residues"
#		print path_file_asa, "path file asa"
#	
#	for lr in liste_residus_pocket:
#		for atom_asa in asa_parsed : 
#			if atom_asa["resSeq"] == lr : 
#				# calcul SASA -> accessibility / (VdW atom) ^2 cf (PETIT 2009)
#				value_ac = atom_asa["ABS"]/(4*pi * (1+atom_asa["radius"]) * (1+atom_asa["radius"]) )
#				if value_ac >=0.2 : 
#					list_atom_retrieve.append (atom_asa["atomSeq"])
#	
#	
#	return list_atom_retrieve
				
# calcul de la derive du log ... que la pente utile en fait





	

# recupere la liste intersection de 2 listes
#def intersection(liste1,liste2):
#	inter = []
#	for i1 in liste1:
#		for i2 in liste2:
#			if i1 == i2:
#				if not i1 in inter : 
#					inter.append(i1)
#	return inter

#def ligand_supp(fichier):
#	# atom count
#	os.system("awk 'NF==9 {print $0}' ligand.mol | wc -l > temp_atom.txt")
#	file = open("temp_atom.txt",'r')
#	n = file.readline()
#	n = string.replace(n,'\n','')
#	# une ligne en trop
#	n = int(n)-1
#	file.close()
#	fichier.write("ligand_atom_count\t"+str(n)+"\n")
#	# bond count
#	os.system("awk 'NF==6 {print $0}' ligand.mol | wc -l > temp_bond.txt")
#	file = open("temp_bond.txt",'r')
#	n = file.readline()
#	n = string.replace(n,'\n','')
#	file.close()
#	fichier.write("ligand_bond_count\t"+n+"\n")






######################################## FONCTIONS TEMPORAIRES A UTILISER, MODIFIER, SUPPRIMER #################################################
#def hiden(extend):
#	# accessibilite a 5 angs pour le gap
#	os.system("naccess -p 2.5 protein.pdb > /dev/null")
#	os.system("mv protein.rsa naccess/protein_5.rsa")
#	os.system("rm protein.asa")
#	os.system("rm protein.log")
#	# nombre de residus accessibles a 1 angs
#	os.system("more protomol/acc_residues_around_"+str.replace(str(extend),'.','_')+".rsa | wc -l > temp_1.out")
#	file = open("temp_1.out",'r')
#	nb_1 = file.readline()
#	nb_1 = int(string.replace(nb_1,'\n',''))
#	file.close()
#	os.system("rm temp_1.out")
#	#print nb_1
#	# les residus a x angs du protomol
#	linesRES = PDB("protomol/residues_around_"+str.replace(str(extend),'.','_')+".pdb")
#	# liste contenant les numeros des residus a x angs du protomol
#	numRES = []
#	for lineRES in linesRES:
#		numRES.append(int(lineRES.rNum()))
#	# on trie la liste pour plus de facilite pour la suite
#	numRES.sort()
#	#num = intersection(numRES,numACC)
#	# le fichier .rsa a 5 angs
#	FileRSA = open("naccess/protein_5.rsa",'r')
#	linesRSA = FileRSA.readlines()
#	FileRSA.close()
#	# on fait un fichier contenant l'accessibilite a 5 angs des residus accessibles de la poche
#	FileOUT = open("naccess/acc_residues_around_"+str.replace(str(extend),'.','_')+"_protein_5.rsa",'w')
#	for n in numRES:
#		for lineRSA in linesRSA:
#			words = string.split(lineRSA)
#			if words[0] == "RES":
#				if len(words)==14:
#					if int(words[3])==n and float(words[5])>=20.:
#						FileOUT.write(str(lineRSA))
#						break
#				elif len(words)==13:
#					if int(words[2])==n and float(words[4])>=20.:
#						FileOUT.write(str(lineRSA))
#						break
#				else:
#					print "probleme fichier naccess/acc_residues.rsa"
#	FileOUT.close()
#	# nombre de residus accessibles a 5 angs
#	os.system("more naccess/acc_residues_around_"+str.replace(str(extend),'.','_')+"_protein_5.rsa | wc -l > temp_5.out")
#	file = open("temp_5.out",'r')
#	nb_5 = file.readline()
#	nb_5 = int(string.replace(nb_5,'\n',''))
#	file.close()
#	os.system("rm temp_5.out")
#	print nb_5
#	print nb_5/float(nb_1)

#def buried_function(ray):
#	#ray = abs(ray)
#	Ro = 2.
#	Do = 1.
#	if ray <= 2.:
#		return 1.
#	return exp(-((ray-Ro)**2)/(Do**2))

#def get_8_angs():
#	memory = []
#	chain = []
#	extend = 8.
#	FileMOL = open("protomol/protomol.pdb",'r')
#	linesMOL = FileMOL.readlines()
#	FileMOL.close()
#	prot = PDB("protein.pdb")
#	for lineMOL in linesMOL:
#		wordsMOL = string.split(lineMOL)
#		if wordsMOL[0] == "HETATM":
#			xMOL = float(wordsMOL[5])
#			yMOL = float(wordsMOL[6])
#			zMOL = float(wordsMOL[7])
#			for linesPDB in prot:
#				if linesPDB.rType() == "AMINO-ACID":
#					for linePDB in linesPDB:
#						wordsPDB = string.split(linePDB.crds())
#						xPDB = float(wordsPDB[0])
#						yPDB = float(wordsPDB[1])
#						zPDB = float(wordsPDB[2])
#						dist = sqrt((xMOL-xPDB)**2+(yMOL-yPDB)**2+(zMOL-zPDB)**2)
#						if dist <= extend:
#							flag = 0
#							for m in memory:
#								if int(linesPDB.rNum()) == m:
#									flag = 1
#									break
#							if flag == 0:
#								memory.append(int(linesPDB.rNum()))
#								chain.append(linesPDB.chnLbl())
#							break
#	#memory.sort()	# pour pouvoir trier il faudrait faire un dico ... car ici la chaine correspond plus au res
#	FileOUT = open("protomol/residues_around_"+str.replace(str(extend),'.','_')+".pdb",'w')
#	if(len(memory)==0):
#		print "poche trop loin de la proteine => a supprimer du jeu"
#	for i in range(len(memory)):
#		for linesPDB in prot:
#			# si on prend le residu en entier ... sinon ben faut faire autrement
#			if(int(linesPDB.rNum())==memory[i] and linesPDB.chnLbl()==chain[i]):
#				FileOUT.write(str(linesPDB))
#	FileOUT.close()
#
#def hiden_bury(pdb):
#	#get_8_angs()
#	os.system("cp protomol/residues_around_16_0.pdb protomol/residues_around_16_0_protomol_merged.pdb")
#	os.system("more protomol/protomol.pdb >> protomol/residues_around_16_0_protomol_merged.pdb")
#	#lines1 = PDB("protomol/protomol_modif.pdb")
#	os.system("awk 'NF==12 && $12!=\"H\" || NF==11 && $11!=\"H\" {print $0}' protomol/residues_around_16_0_protomol_merged.pdb > protomol/residues_around_16_0_protomol_merged_modif.pdb")
#	lines1 = PDB("protomol/residues_around_16_0_protomol_merged_modif.pdb")
#	dico_sonde = []
#	n = 0
#	# on recupere tout et on renumerote
#	for res1 in lines1:
#		for at1 in res1:
#			sonde = {}
#			sonde['num'] = n
#			sonde['header'] = at1.header()
#			crds = string.split(at1.crds())
#			sonde['crdX'] = float(crds[0])
#			sonde['crdY'] = float(crds[1])
#			sonde['crdZ'] = float(crds[2])
#			n += 1
#			dico_sonde.append(sonde)
#	#print dico_sonde
#	l = len(dico_sonde)
#	#print l
#	# on calcul les bc
#	for i in range(l):
#		# liste qui va contenir les num des sondes a moins de 4 angs
#		liste = []
#		# liste qui va contenir les distances
#		dist = []
#		# bc
#		bc = 0
#		xi = dico_sonde[i]['crdX']
#		yi = dico_sonde[i]['crdY']
#		zi = dico_sonde[i]['crdZ']
#		for j in range(l):
#			if j != i:
#				xj = dico_sonde[j]['crdX']
#				yj = dico_sonde[j]['crdY']
#				zj = dico_sonde[j]['crdZ']
#				d = sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
#				if d <= 8.:
#					liste.append(j)
#					dist.append(d)
#					bc += 1
#		dico_sonde[i]['liste'] = liste
#		dico_sonde[i]['dist'] = dist
#		dico_sonde[i]['bc'] = bc
#	# on calcule les pw que sur les sondes (mais donc en tenant compte de la proteine)
#	liste_pw = []
#	for i in range(l):
#		if dico_sonde[i]['header'] == 'HETATM':
#			pw = 0
#			n = 0
#			for j in dico_sonde[i]['liste']:
#				pw += dico_sonde[j]['bc']*buried_function(dico_sonde[i]['dist'][n])
#				n += 1
#			liste_pw.append(pw)
#	#print dico_sonde[60]
#	monfile = open("/home/mti-priv91/perot/These/Projet1_Biochemo/Resultats/Bury/distribution.out",'a')
#	monfile.write(pdb+"\n")
#	m = min(liste_pw)
#	#monfile.write("min :"+str(m)+"\n")
#	M = max(liste_pw)
#	#monfile.write("max :"+str(M)+"\n")
#	av = mean(liste_pw)
#	#monfile.write("moyenne :"+str(av)+"\n")
#	#monfile.write("len :"+str(len(liste_pw))+"\n")
#	#monfile.write("max / n :"+str(M/len(liste_pw))+"\n")
#	monfile.write("d."+pdb+" = c(")
#	for i in liste_pw:
#		monfile.write(str("%2.f"%i)+",")
#	monfile.write(")\n")
#	monfile.close()
#
#
#def hiden_bury_1():
#	# tire de PASS (Brady et Stouten)
#	os.system("awk '$3==\"C\" || $3==\"N\" || $3==\"O\" {print $0}' protomol/protomol.pdb > protomol/protomol_modif.pdb")
#	lines1 = PDB("protomol/protomol_modif.pdb")
#	#lines1 = PDB("protomol/protein_protomol_merged.pdb")
#	lines2 = PDB("protomol/protomol_modif.pdb")
#	#lines2 = PDB("protomol/protein_protomol_merged.pdb")
#	#liste_bc = []
#	liste_pw = []
#	dico_sonde = []
#	for res1 in lines1:
#		for at1 in res1:
#			bc = 0
#			sonde = {}
#			liste = []
#			dist = []
#			coord1 = string.split(at1.crds())
#			x1 = float(coord1[0])
#			y1 = float(coord1[1])
#			z1 = float(coord1[2])
#			k = 0
#			for res2 in lines2:
#				for at2 in res2:
#					k += 1
#					coord2 = string.split(at2.crds())
#					x2 = float(coord2[0])
#					y2 = float(coord2[1])
#					z2 = float(coord2[2])
#					d = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
#					if d <= 8.:
#						if res1.rNum() != res2.rNum():
#							bc += 1
#							# moins 1 car apres on commence a zero
#							liste.append(int(res2.rNum())-1)
#							#liste.append(k)
#							dist.append(d)
#			sonde['bc'] = bc
#			sonde['liste'] = liste
#			sonde['dist'] = dist
#			sonde['header'] = at1.header()
#			#liste_bc.append(bc)
#			dico_sonde.append(sonde)
#	for i in range(len(dico_sonde)):
#		if dico_sonde[i]['header'] == "HETATM":
#			pw = 0
#			n = 0
#			for j in dico_sonde[i]['liste']:
#				pw += dico_sonde[j]['bc']*buried_function(dico_sonde[i]['dist'][n])
#				n += 1
#			liste_pw.append(pw)
#	m = min(liste_pw)
#	#M = max(liste_pw)
#	#print M
#
#def graphe():
#	pdb = PDB("protomol/acc_residues_around_4_0.pdb")
#	dico = []
#	for res in pdb:
#		residu = {}
#		for at in res:
#			if at.atmName()=='CA':
#				coord = at.crds()
#				residu['crds'] = coord
#				residu['name'] = at.resName()
#				#print residu
#		dico.append(residu)
#	#print dico
#	for i in range(len(dico)-1):
#		crdsi = dico[i]['crds']
#		ci = string.split(crdsi)
#		xi = float(ci[0])
#		yi = float(ci[1])
#		zi = float(ci[2])
#		for j in range((i+1),len(dico)):
#			crdsj = dico[j]['crds']
#			cj = string.split(crdsj)
#			xj = float(cj[0])
#			yj = float(cj[1])
#			zj = float(cj[2])
#			d = sqrt((xi-xj)**2+(yi-yj)**2+(zi-zj)**2)
#			if d <= 5.:
#				print dico[i]['name'],dico[j]['name'],d
#
#
#def info1_ligand(fileout,info,lines):
#	fileout.write("ligand_"+string.lower(info)+"\t")
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)>1 and words[0] == info:
#			if info=="POTENCY":
#				if words[1] != "NA":
#					words[1] = str("%.2f" %log(float(words[1])))
#			fileout.write(words[1]+"\n")
#			na = 0
#			break
#	if na==1:
#		fileout.write("NA\n")
#
#
#def info2_ligand(fileout,info1,info2,lines):
#	fileout.write("ligand_"+string.lower(info1)+"_"+string.lower(info2)+"\t")
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)>2 and words[0] == info1 and words[1] == info2:
#			fileout.write(words[2]+"\n")
#			na = 0
#			break
#	if na==1:
#		fileout.write("NA\n")
#
#
#def info3_ligand(fileout,lines):
#	fileout.write("ligand_chiral_atom_count\t")
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)>3 and words[0] == "Chiral" and words[1] == "atom":
#			fileout.write(words[3]+"\n")
#			na = 0
#			break
#	if na==1:
#		fileout.write("NA\n")
#	fileout.write("ligand_aromatic_bond_count\t")
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)==4 and words[0] == "Aromatic" and words[1] == "bond" and words[2] == "count":
#			fileout.write(words[3]+"\n")
#			na = 0
#			break
#	if na==1:
#		fileout.write("NA\n")
#
#
#def energy(fileout,info,lines):
#	fileout.write("ligand_"+string.lower(info)+"_energy"+"\t")
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)>=4 and words[0] == info and words[1] == "Ligand" and words[2] == "interactions":
#			fileout.write(words[3]+"\n")
#			na = 0
#			break
#	if na==1:
#		fileout.write("NA\n")
#
#
#def ligand(fileout):
#	info = open("info/info.txt",'r')
#	lines_info = info.readlines()
#	info.close()
#	info1_ligand(fileout,"miLogP",lines_info)
#	info1_ligand(fileout,"TPSA",lines_info)
#	# nombre d'atomes lourds
#	info1_ligand(fileout,"natoms",lines_info)
#	info1_ligand(fileout,"MW",lines_info)
#	info1_ligand(fileout,"nON",lines_info)
#	info1_ligand(fileout,"nOHNH",lines_info)
#	info1_ligand(fileout,"nviolations",lines_info)
#	info1_ligand(fileout,"nrotb",lines_info)
#	# voir FaFDrug ... RPBS ... en attente
#	#info2_ligand(fileout,"Formal","charge",lines_info)
#	# nombre d'atomes (lourds et hydrogenes)
#	info2_ligand(fileout,"Atom","count",lines_info)
#	info2_ligand(fileout,"Bond","count",lines_info)
#	info3_ligand(fileout,lines_info)
#	info1_ligand(fileout,"POTENCY",lines_info)
#	energy(fileout,"Internal",lines_info)
#
#
#def info_energy(fileout,info1,info2,lines):
#	fileout.write("interaction_"+info1+"_"+info2+"\t")
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)>=3 and words[1] == info1 and words[2] == info2:
#			if words[0] == "no" or words[0] == "zero":
#				fileout.write("0\n")
#			elif words[0] == "yes" or words[0] == "one" or words[0] == "1" or words[0] == "some":
#				fileout.write("1\n")
#			else:
#				fileout.write(words[0]+"\n")
#			na = 0
#			break
#	# si rien alors c'est 0
#	if na==1:
#		fileout.write("0\n")

#
#def tension(fileout,lines):
#	fileout.write("interaction_ligand_tension\t")
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)>3 and words[0] == "100" and (string.lower(words[1]) == "steps" or string.lower(words[1]) == "step"):
#			step100 = float(words[2])
#			#fileout.write(words[2]+"\n")
#			na = 0
#			break
#		elif len(words)>3 and words[0] == "100":
#			step100 = float(words[1])
#			#fileout.write(words[1]+"\n")
#			na = 0
#			break
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)>=3 and words[0] == "random" and words[1] == "search":
#			randomsearch = float(words[2])
#			tension = "%.2f" %(step100 - randomsearch)
#			fileout.write(str(tension)+"\n")
#			na = 0
#			break
#		elif len(words)>=3 and words[0] == "lowest" and words[1] == "energy":
#			randomsearch = float(words[4])
#			tension = "%.2f" %(step100 - randomsearch)
#			fileout.write(str(tension)+"\n")
#			na = 0
#			break
#	if na==1:
#		fileout.write("NA\n")
#
#
#def energy_interactions(fileout,info,lines):
#	fileout.write("interaction_ligand_"+string.lower(info)+"\t")
#	na = 1
#	for l_info in lines:
#		words = string.split(l_info)
#		if len(words)>=5 and words[0] == info and words[1] == "-" and words[2] == "Ligand" and words[3] == "interactions":
#			fileout.write(words[4]+"\n")
#			na = 0
#			break
#	if na==1:
#		fileout.write("NA\n")
#
#def interaction(fileout):
#	info = open("info/info.txt",'r')
#	lines_info = info.readlines()
#	info.close()
#	info_energy(fileout,"cation","pi",lines_info)
#	info_energy(fileout,"salt","bridge",lines_info)
#	info_energy(fileout,"aromatic","stacking",lines_info)
#	tension(fileout,lines_info)
#	energy_interactions(fileout,"Protein",lines_info)
#	energy_interactions(fileout,"Water",lines_info)
#
#





# FONCTION : CA(x) et SA(x) et RA(x) et PLB voir Soga et al.
#def soga(seqA,seqB,fileout):
#def soga(fileout,pdb):
#	if pdb != "1hwi" and pdb != "1n1m":
#		#nameAcc = "protomol/acc_residues_around_"+str.replace(str(extend),'.','_')
#		#pocketAcc = PDB(nameAcc+".pdb")
#		proteinAcc = open("naccess/acc_residues.rsa",'r')
#		lineAcc = proteinAcc.readlines()
#		proteinAcc.close()
#		numAcc = []
#		for line in lineAcc:
#			words = string.split(line)
#			numAcc.append(int(words[3]))
#		proteinAccPDB = open("naccess/acc_residues.pdb",'w')
#		prot = PDB("protein.pdb")
#		for resprot in prot:
#			for n in numAcc:
#				if n == int(resprot.rNum()):
#					proteinAccPDB.write(str(resprot))
#					break
#		proteinAccPDB.close()
#		nameAcc = "naccess/acc_residues"
#		protAcc = PDB(nameAcc+".pdb")
#		nameRes = "protomol/residues_around_"+str.replace(str(extend),'.','_')
#		pocketRes = PDB(nameRes+".pdb")
#		seqB = protAcc.aaseq()
#		seqA = pocketRes.aaseq()
#		sA = len(seqA)
#		sB = len(seqB)
#		PLB = 0.
#		for aa in AAun:
#			nbA = string.count(seqA,aa)
#			CAx = nbA/float(sA)
#			nbB = string.count(seqB,aa)
#			SAx = nbB/float(sB)
#			if SAx == 0.:
#				RAx = 0.
#			else:
#				RAx = CAx/SAx
#			PLB += nbA*RAx
#			fileout.write(str(SAx)+" ")
#			#fileout.write("%.2f " % CAx+" %.2f " % SAx+" %.2f " % RAx+" ")
#			#fileout.write("pocket_CA_"+aa+"\t%.2f " % CAx+"\n")
#			#fileout.write("pocket_SA_"+aa+"\t%.2f " % SAx+"\n")
#			#fileout.write("pocket_RA_"+aa+"\t%.2f " % RAx+"\n")
#		#fileout.write("pocket_PLB\t%.2f " % PLB+"\n")
#		fileout.write("\n")

#def resolution(fileout):
#	code = string.upper(pdb)
#	code = code+".pdb"
#	os.system("more "+code+" | grep \"RESOLUTION\.\" > temp_res.out")
#	file = open("temp_res.out",'r')
#	res = file.readline()
#	res = string.replace(res,'\n','')
#	file.close()
#	resolution = string.split(res)[3]
#	fileout.write(resolution+"\n")
#
#def olivier(pdb):
#	os.system("mkdir /home/mti-priv91/perot/Desktop/Ligands/"+pdb)
#	os.system("cp ligand.mol /home/mti-priv91/perot/Desktop/Ligands/"+pdb)
#	os.system("cp ligand.pdb /home/mti-priv91/perot/Desktop/Ligands/"+pdb)



# FONCTION : calcul du volume de la surface et du rapport ie compacite a partir de msms
#def volume_msms(flag,path,fileout):
#	os.chdir("protomol")
#	# recuperation du fichier atmtypenumbers ... necessaire pour pdb_to_xyzr
#	os.system("cp -f "+path+"/These/Projet1_Biochemo/Codes/Python/atmtypenumbers .")
#	# creation du fichier .xyzr
#	os.system("pdb_to_xyzr protomol.pdb > protomol.xyzr")
#	# calcul du volume et de la surface dans le fichier .msms
#	os.system("msms -if protomol.xyzr > protomol.msms")
#	# recuperation du volume
#	os.system("awk '$1==\"Total\" && $2==\"ses_volume:\" {print $3}' protomol.msms > temp_vol.out")
#	fileV = open("temp_vol.out",'r')
#	nV = fileV.readline()
#	nV = string.replace(nV,'\n','')
#	fileV.close()
#	os.system("rm temp_vol.out")
#	# recuperation de la surface a partir de la valeur du volume
#	os.system("awk '$3==\""+nV+"\" && NF==4 {print $4}' protomol.msms > temp_surf.out")
#	fileS = open("temp_surf.out",'r')
#	nS = fileS.readline()
#	nS = string.replace(nS,'\n','')
#	fileS.close()
#	os.system("rm temp_surf.out")
#	nV = float(nV)
#	nS = float(nS)
#	# compacite
#	compacity = nV/nS
#	compacity = "%.2f" %compacity
#	nV = int(nV)
#	nS = int(nS)
#	os.chdir("..")
#	fileout.write(flag+"_volume\t"+str(nV)+"\n")
#	fileout.write(flag+"_surface\t"+str(nS)+"\n")
#	fileout.write(flag+"_compacity\t"+str(compacity)+"\n")


# FONCTION : calcul du volume de la surface et du rapport ie compacite a partir de chimera
#def compute_chimera_hide(flag,path,proto):
#	if flag=="pocket":
#		# cp du protomol pour avoir le meme nom pour chimera
#		os.system("cp protomol/"+proto+"-protomol.mol2 protomol/protomol_chimera_temp.mol2")
#		#print "awk '$2==\"C\" {print $0} $2==\"N\" {print $0}' protomol_chimera_temp.mol2 > protomol_chimera.mol2"
#		#os.system("awk '$2==\"C\" {print $0} $2==\"N\" {print $0} $2==\"O\"' protomol/protomol_chimera_temp.mol2 > protomol/protomol_chimera.mol2")
#		file1 = open("protomol/protomol_chimera_temp.mol2",'r')
#		file2 = open("protomol/protomol_chimera.mol2",'w')
#		lines1 = file1.readlines()
#		file1.close()
#		for line1 in lines1:
#			words = string.split(line1)
#			if len(words) == 0 or len(words) == 1:
#				file2.write(line1)
#			else:
#				if words[1] != "H":
#					file2.write(line1)
#		file2.close()
#	# flag p pour poche, l pour ligand
#	# redirection du message de chimera pour la registration dans le /dev/null
#	#os.system("chimera --nogui < "+path+"/These/Projet1_Biochemo/Codes/Python/get_chimera_"+flag+".txt > descriptors/chimera_volume_"+flag+".out 2> /dev/null")
#	os.system("chimera --nogui < "+path+"/Pocket_project/src/get_chimera_"+flag+".txt > descriptors/chimera_volume_"+flag+".out 2> /dev/null") ##LR
#	os.system("awk '$1==\"MSMS\" && $6==\"volume\" {print $8}' descriptors/chimera_volume_"+flag+".out > temp_vol.out")
#	fileV = open("temp_vol.out",'r')
#	nV = fileV.readline()
#	nV = string.replace(nV,'>\n','')
#	fileV.close()
#	os.system("rm temp_vol.out")
#	os.system("awk '$1==\"MSMS\" && $6==\"area\" {print $8}' descriptors/chimera_volume_"+flag+".out > temp_surf.out")
#	fileS = open("temp_surf.out",'r')
#	nS = fileS.readline()
#	nS = string.replace(nS,'>\n','')
#	fileS.close()
#	os.system("rm temp_surf.out")
#	nV = float(nV)
#	nS = float(nS)
#	# compacite
#	compacity = nV/nS
#	compacity = "%.2f" %compacity
#	nV = int(nV)
#	nS = int(nS)

#def volume_chimera_hide(path):
#	compute_chimera_hide("pocket",path)
#	if (os.path.isfile("ligand.mol")!=0):
#		compute_chimera_hide("ligand",path)
#
#def volume_msms_hide(path):
#	os.chdir("protomol")
#	os.system("awk '$3==\"C\" {print $0} $3==\"N\" {print $0}' protomol.pdb > protomol_new.pdb")
#	# recuperation du fichier atmtypenumbers ... necessaire pour pdb_to_xyzr
#	os.system("cp -f "+path+"/These/Projet1_Biochemo/Codes/Python/atmtypenumbers .")
#	# creation du fichier .xyzr
#	os.system("pdb_to_xyzr protomol_new.pdb > protomol_new.xyzr")
#	# calcul du volume et de la surface dans le fichier .msms
#	os.system("msms -if protomol_new.xyzr > protomol.msms")
#	# recuperation du volume
#	os.system("awk '$1==\"Total\" && $2==\"ses_volume:\" {print $3}' protomol.msms > temp_vol.out")
#	fileV = open("temp_vol.out",'r')
#	nV = fileV.readline()
#	nV = string.replace(nV,'\n','')
#	fileV.close()
#	os.system("rm temp_vol.out")
#	# recuperation de la surface a partir de la valeur du volume
#	os.system("awk '$3==\""+nV+"\" && NF==4 {print $4}' protomol.msms > temp_surf.out")
#	fileS = open("temp_surf.out",'r')
#	nS = fileS.readline()
#	nS = string.replace(nS,'\n','')
#	fileS.close()
#	os.system("rm temp_surf.out")
#	nV = float(nV)
#	nS = float(nS)
#	# compacite
#	compacity = nV/nS
#	compacity = "%.2f" %compacity
#	nV = int(nV)
#	nS = int(nS)
#	os.chdir("..")


#def atom_ligand(fichier,fileout,na):
#	os.system("grep \"HETATM\" ligand"+fichier+".pdb | wc -l > temp_ligand_1.out")
#	fileA = open("temp_ligand_1.out",'r')
#	nA = fileA.readline()
#	nA = string.replace(nA,'\n','')
#	fileA.close()
#	os.system("rm temp_ligand_1.out")
#	fileout.write("ligand_atom_counts\t"+nA+"\n")
#
#
#def bond_ligand(fichier,fileout,na):
#	os.system("grep \"CONECT\" ligand"+fichier+".pdb | wc -l > temp_ligand.out")
#	fileB = open("temp_ligand.out",'r')
#	nB = fileB.readline()
#	nB = string.replace(nB,'\n','')
#	fileB.close()
#	os.system("rm temp_ligand.out")
#	fileout.write("ligand_bond_counts\t"+nB+"\n")



#def atomes_codes(flag,name,fileout):
#	if flag == "pocket":
#		#print name
#		filePRO = os.listdir("protomol/")
#		nom = string.split(name,'/')[-1]
#		for fi in filePRO:
#			# on efface le fichier preexistant
#			if re.search(nom+".mol2",fi):
#				os.system("rm protomol/"+fi)
#		pdb = PDB(name+".pdb")
#		liste_res = []
#		for res in pdb:
#			liste_res.append(res.rNum())
#		#print liste_res
#		for i in liste_res:
#			#print i
#			#print "awk '$7==\""+i+"\" {print $0}' protein.mol2"
#			os.system("awk '$7==\""+i+"\" {print $0}' protein.mol2 >> "+name+".mol2")
#	n_tot = 0
#	# voir article Lu et al. J. Chem. Inf. Model. 2007
#	for at in ["C.3","C.2","C.ar","C.cat","N.3","N.2","N.ar","N.am","N.4","N.pl3",
#	"O.3","O.2","O.co2","S.3","S.2","S.O","S.O2","P.3","F","Cl","H"]:
#		# atomes du fichier
#		os.system("awk '$6==\""+at+"\"' "+name+".mol2 | wc -l > temp.out")
#		file = open("temp.out",'r')
#		n = file.readline()
#		n = string.replace(n,'\n','')
#		n_tot += int(n)
#		file.close()
#		os.system("rm temp.out")
#		fileout.write(flag+"_"+at+"\t"+n+"\n")
#	fileout.write(flag+"_atomes\t"+str(n_tot)+"\n")

# pas interessant ici car en fait je veux la distance relative donc avec le signe
# Plane.distanceFrom me donne la distance en terme de distance donc en valeur absolue !!
#def test_planarity(filein,fileout):
#	params = least_square_plane(filein)
#	A = params[0]
#	B = params[1]
#	C = params[2]
#	p1 = Vector(0.,0.,C)
#	p2 = Vector(0.,1.,B+C)
#	p3 = Vector(1.,0.,A+C)
#	n = Vector(A,B,-1)
#	plan = Plane(p1,p2,p3)
#	dist = []
#	pdb = PDB(filein)
#	for res1 in pdb:
#		for at1 in res1:
#			xyz1 = at1.xyz()
#			x1 = xyz1[0]
#			y1 = xyz1[1]
#			z1 = xyz1[2]
#			point = Vector(x1,y1,z1)
#			d = Plane.distanceFrom(plan,point)
#			dist.append(d)
#	d1 = abs(min(dist))
#	d2 = max(dist)
#	d1d2 = d1+d2
#	d1d2 = "%.2f" %d1d2
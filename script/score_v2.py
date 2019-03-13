import sys
#pathsrc= "/home/hecini/Research/stage_HECINI/script/"   #direcrtory ou vous avez les script PDB6.py 
sys.path.append('/home/hecini/Research/stage_HECINI/script/')
from PDB6 import *

import string
from glob import glob # pour lrie l ensemble des fichier en un seul coup 
import pandas as pd # pour transformer le dictionnaire en dataFrame






#lecture des fichiers et creation d un dictionnaire 

pdb_files = glob('/home/hecini/Research/stage_HECINI/data/PR2/pockets/*.pdb')
liste_finale = []
dic = {}

for fileName in pdb_files:
	structure_id = fileName.rsplit('/', 1)[1][:-5]
	pdb_obj = PDB(fileName)
	liste = []
	for res in pdb_obj : 
	   for atm in res:
	       atmNm = atm.atmName()
	       resNum = atm.resNum()
	       ch = atm.chnLbl()
	       atom = string.join([resNum,ch,atmNm],"_")
	       liste.append(atom)
	liste_finale.append(liste)
	dic[structure_id] = liste




#### calcul des scores : 

dic_score = {}

for k in range(0,len(dic.keys())):
	for l in range(k+1,len(dic.keys())):
		nbrIntersect = len(set(dic.values()[k]).intersection(dic.values()[l]))
		nbrAtom1 = len(dic.values()[k])
		nbrAtom2 = len(dic.values()[l])
		y = float(nbrIntersect)/(nbrAtom1+nbrAtom2-nbrIntersect)
		Sc = float("{0:.4f}".format(y))
		dic_score[dic.keys()[k]+';'+dic.keys()[l]]= Sc



### transformation en data frame --> csv 

dt = pd.DataFrame (dic_score.items(), columns=['poches', 'score'])

dt.to_csv('score3.csv')





#### avec liste de liste : 

"""liste_score = []

for i in range(0,len(liste_finale)-1):
	for j in range(i+1,len(liste_finale)):
		n_inter = len(set(liste_finale[i]).intersection(liste_finale[j]))
		nbrAtom3 = len(liste_finale[i])
		nbrAtom4 = len(liste_finale[j])
		x2 = float(n_inter)/(nbrAtom3+nbrAtom4-n_inter)
		SO2 = float("{0:.4f}".format(x2))
		liste_score.append(SO2)"""

#### calcul des scores : 

"""dic_score = {}

for poche in dic:
	for poche2 in dic:
		nbrIntersect = len(set(dic[poche]).intersection(dic[poche2]))
		nbrAtom1 = len(dic[poche])
		nbrAtom2 = len(dic[poche2])
		x = float(nbrIntersect)/(nbrAtom1+nbrAtom2-nbrIntersect)
		SO = float("{0:.4f}".format(x))
		dic_score[poche+'_'+poche2]= SO"""











     




import sys
#pathsrc= "/home/hecini/Research/stage_HECINI/script/"   #direcrtory ou vous avez les script PDB6.py 
sys.path.append('/home/hecini/Research/stage_HECINI/script/')
from PDB6 import *
import string
from glob import glob

pdb_files = glob('/home/hecini/Research/stage_HECINI/data/PR2/pockets/*.pdb')


liste_finale = []
dic = {}

for fileName in pdb_files:
	structure_id = fileName.rsplit('/', 1)[1][0:13]
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

dic_score = {}

for poche in dic:
	for poche2 in dic:
		nbrIntersect = len(set(dic[poche]).intersection(dic[poche2]))
		nbrAtom1 = len(dic[poche])
		nbrAtom2 = len(dic[poche])
		x = float(nbrIntersect)/(nbrAtom1+nbrAtom2-nbrIntersect)
		SO = float("{0:.4f}".format(x))
		dic_score[poche+poche2]= SO



import pandas as pd 

dt = pd.DataFrame (dic_score.items(), columns=['poches', 'score'])

dt.to_excel('output1.xlsx', engine='xlsxwriter') 


     


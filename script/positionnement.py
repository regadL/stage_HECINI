import sys
sys.path.append('/home/hecini/Research/stage_HECINI/script/')
from PDB6 import *

import string
from glob import glob # pour lrie l ensemble des fichier en un seul coup 
import pandas as pd # pour transformer le dictionnaire en dataFrame
import operator 



#lecture des fichiers et creation d un dictionnaire 

pdb_files = glob('/home/hecini/Research/stage_HECINI/data/PR2/positionnement3D_PDB/*.pdb')

liste = []
dic = {}

for fileName in pdb_files:
	structure_id = fileName.rsplit('/', 1)[1][:] # recuperer le nom de chaque poche 
	pdb_obj = PDB(fileName) # le contenu de chaque poche
	somme = (0,0,0)
	occurence = 0
	for pdb in pdb_obj : #lecture ligne par ligne 
		for atm in pdb:
			#if (atm.atmName() == 'CA'):
			occurence= occurence+1
			somme = tuple(map(operator.add, somme,atm.xyz()))
	liste.append(somme)
	barycentre = tuple(cordonne/occurence for cordonne in somme)
	
	dic[structure_id] = list(barycentre)



### transformation en data frame --> csv 

dt = pd.DataFrame (dic.items(), columns=['poche', 'barycentre'])
#df3 = pd.DataFrame(dt['barycentre'].values.tolist(), columns=['x','y','z'])

dt[['x','y','z']] = pd.DataFrame(dt.barycentre.values.tolist(), index= dt.index)
dt.to_csv('barycentre_new')








#liste_finale.append(liste)
#dic[structure_id] = liste
#print(len(liste))
#print(liste[0])

#print(dic.keys()[1])
#print( dic.values()[1])

#print(liste_finale[1])




	   

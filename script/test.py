import sys
sys.path.append('/home/hecini/Research/stage_HECINI/script/')
from PDB6 import *

import string
from glob import glob # pour lrie l ensemble des fichier en un seul coup 
import pandas as pd # pour transformer le dictionnaire en dataFrame
import operator 



#lecture des fichiers et creation d un dictionnaire 
liste= []
dic = {}

pdb_obj = PDB("/home/hecini/Research/stage_HECINI/data/PR2/pockets/1HSI_pocket0_atm.pdb") # le contenu de chaque poche
"""liste = []
somme = (0,0,0)
for pdb in pdb_obj :   
	#print(liste)
	#dic[atm.atmName()] = somme
	somme = (0,0,0)
	for atm in pdb:
		if (atm.atmName() == 'CA'):
			somme = tuple(map(operator.add, somme,atm.xyz()))
			#atom = string.join([atm.atmName(),cord],"_"
	liste.append(somme)"""

somme = (0,0,0)

for pdb in pdb_obj : #lecture ligne par ligne 
	#dic[atm.atmName()] = somme
	for atm in pdb:
		if (atm.atmName() == 'CA'):
			print(atm.xyz())
			somme = tuple(map(operator.add, somme,atm.xyz()))
liste.append(somme)

print(count('CA'))
print(liste)
	   

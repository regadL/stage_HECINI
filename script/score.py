

import sys
#pathsrc= "/home/hecini/Research/stage_HECINI/script/"   #direcrtory ou vous avez les script PDB6.py 
sys.path.append('/home/hecini/Research/stage_HECINI/script/')
from PDB6 import *
import string
from glob import glob

pdb_files = glob('/home/hecini/Research/stage_HECINI/data/PR2/pockets/*.pdb')
print(pdb_files)

#for fileName in pdb_files:
	
     



pdb_obj = PDB("/home/hecini/Research/stage_HECINI/data/PR2/pockets/1IDB_pocket0_atm.pdb")

list1 = []

for res in pdb_obj : 
   for atm in res:
       atmNm = atm.atmName()
       resNum = atm.resNum()
       ch = atm.chnLbl()
       atom = string.join([resNum,ch,atmNm],"_")
       list1.append(atom)



       

pdb_obj1 = PDB("/home/hecini/Research/stage_HECINI/data/PR2/pockets/1HSI_pocket0_atm.pdb")
list2 = []

for res in pdb_obj1 : 
   for atm in res:
       atmNm = atm.atmName()
       resNum = atm.resNum()
       ch = atm.chnLbl()
       atom = string.join([resNum,ch,atmNm],"_")
       list2.append(atom)




       

nbrIntersect = len(set(list2).intersection(list1))
nbrAtom1 = len(list1)
nbrAtom2 = len(list2)
SO = float(nbrIntersect)/(nbrAtom1+nbrAtom2-nbrIntersect)
print(SO)

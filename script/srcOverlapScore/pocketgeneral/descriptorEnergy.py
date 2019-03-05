"""
BORREL Alexandre
04-2012
calcul descriptor
"""

import runOtherProg
import string
import os
import parseNACCESS



def ratioOriginal(fichierin, flag, filout, na = 0):
    if na==0:
        point = 1
        if flag == "protomol":
            fichier = fichierin
        elif flag == "ligand":
            fichier = "ligand"+fichierin+".mol"
            point = 0
        elif flag == "pocket":
            fichier = runOtherProg.babelPDBtoMOL2(fichierin)
            
        elif flag == "protein":
            fichier = "protein.mol2"
        else:
            print "flag inconnu"
        if point == 1:
            os.system("grep \" O\.\" "+fichier+" | wc -l > temp_O.out")
            fileO = open("temp_O.out",'r')
            nO = fileO.readline()
            nO = int(string.replace(nO,'\n',''))
            fileO.close()
            os.system("rm temp_O.out")
            os.system("grep \" N\.\" "+fichier+" | wc -l > temp_N.out")
            fileN = open("temp_N.out",'r')
            nN = fileN.readline()
            nN = int(string.replace(nN,'\n',''))
            fileN.close()
            os.system("rm temp_N.out")
            os.system("grep \" S\.\" "+fichier+" | wc -l > temp_S_all.out")
            # sans les S des methionine ... cf Olivier
            os.system("grep \" S\..*CYS\" "+fichier+" | wc -l > temp_S.out")
            fileS = open("temp_S.out",'r')
            nS = fileS.readline()
            nS = int(string.replace(nS,'\n',''))
            fileS.close()
            os.system("rm temp_S.out")
            fileSall = open("temp_S_all.out",'r')
            nSall = fileSall.readline()
            nSall = int(string.replace(nSall,'\n',''))
            fileSall.close()
            os.system("rm temp_S_all.out")
            os.system("grep \" C\.3\" "+fichier+" | wc -l > temp_C.out")
            fileC = open("temp_C.out",'r')
            nC = fileC.readline()
            nC = int(string.replace(nC,'\n',''))
            fileC.close()
            os.system("rm temp_C.out")
        elif point == 0:
            os.system("grep \" O \" "+fichier+" | wc -l > temp_O.out")
            fileO = open("temp_O.out",'r')
            nO = fileO.readline()
            nO = int(string.replace(nO,'\n',''))
            fileO.close()
            os.system("rm temp_O.out")
            os.system("grep \" N \" "+fichier+" | wc -l > temp_N.out")
            fileN = open("temp_N.out",'r')
            nN = fileN.readline()
            nN = int(string.replace(nN,'\n',''))
            fileN.close()
            os.system("rm temp_N.out")
            os.system("grep \" S \" "+fichier+" | wc -l > temp_S.out")
            fileS = open("temp_S.out",'r')
            nS = fileS.readline()
            nS = int(string.replace(nS,'\n',''))
            nSall = nS
            fileS.close()
            os.system("rm temp_S.out")
            os.system("grep \" C \" "+fichier+" | wc -l > temp_C.out")
            fileC = open("temp_C.out",'r')
            nC = fileC.readline()
            nC = int(string.replace(nC,'\n',''))
            fileC.close()
            os.system("rm temp_C.out")
        # polarity ratio -> a revoir plus propre
        try : 
            ratio = (nO+nN+nS)/float(nO+nN+nSall+nC)
            ratio = "%.2f" %ratio
        except :
            ratio = "NA"
        # hydrophobicity ratio voir Burgoyne et Jackson Bioinformatics
        try :
            ratio2 = (nC+nS)/float(nO+nN+nSall+nC)
            ratio2 = "%.2f" %ratio2
        except : 
            ratio2 = "NA"
    elif na==1:
        ratio = "NA"
    filout.write("pocket_polarity_pocket_" + flag + "\t"+ratio+"\n")
    filout.write("pocket_hydrophobicity_pocket_" + flag + "\t"+ratio2+"\n")
    
    return fichier
    

def accessibilityPocketbyAtom (path_file_asa):
    """Retrieve sum of atom ABS
    arg: pah file asa
    return: accessibility (float)"""
    
    list_atom_asa = parseNACCESS.fileASA(path_file_asa)    
    accesibility = 0.0
    
    for atom_asa in list_atom_asa : 
        accesibility = accesibility + atom_asa["ABS"]
    
    return accesibility

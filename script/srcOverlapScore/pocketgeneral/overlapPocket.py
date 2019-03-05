import parsePDB
import parseNACCESS
import superposeStructure
import writePDBfile
import os


def scoreOverlap (path_file_pocket1, path_file_pocket2):
    
    list_atom_pocket1 = parsePDB.loadCoordSectionPDB(path_file_pocket1)
    list_atom_pocket2 = parsePDB.loadCoordSectionPDB(path_file_pocket2)
    
    int_nb_atom = len (list_atom_pocket1) + len (list_atom_pocket2)
    
    
    int_commun = 0.0
    for dico_atom_pocket1 in list_atom_pocket1 : 
        for dico_atom_pocket2 in list_atom_pocket2 : 
            if dico_atom_pocket1["x"] == dico_atom_pocket2["x"] and dico_atom_pocket1["y"] == dico_atom_pocket2["y"] and dico_atom_pocket1["z"] == dico_atom_pocket2["z"] : 
                int_commun  = int_commun + 1
                break
    
    if (int_nb_atom - int_commun) == 0.0 : 
        return 0.0
    else : 
        return int_commun / (int_nb_atom - int_commun)
#     return int_commun


def communAtom(path_file_pocket1, path_file_pocket2) : 
    
    list_atom_pocket1 = parsePDB.loadCoordSectionPDB(path_file_pocket1)
    list_atom_pocket2 = parsePDB.loadCoordSectionPDB(path_file_pocket2)
    
    int_nb_atom = len (list_atom_pocket1) + len (list_atom_pocket2)
    
    
    int_commun = 0.0
    for dico_atom_pocket1 in list_atom_pocket1 : 
        for dico_atom_pocket2 in list_atom_pocket2 : 
            if dico_atom_pocket1["x"] == dico_atom_pocket2["x"] and dico_atom_pocket1["y"] == dico_atom_pocket2["y"] and dico_atom_pocket1["z"] == dico_atom_pocket2["z"] : 
                int_commun  = int_commun + 1
                break
    
    return int_commun , len (list_atom_pocket1), len (list_atom_pocket2)


def MO (path_file_pocket1, path_file_pocket2, path_file_pocket1_asa, path_file_pocket2_asa, debug = 0):
    
    
    if debug : 
        print "***************"
        print path_file_pocket1
        print path_file_pocket2
        print path_file_pocket1_asa
        print path_file_pocket2_asa
        print "***************"
    
    print os.path.split(path_file_pocket1)
    
    l_atom_pocket1 = parsePDB.loadCoordSectionPDB(path_file_pocket1)
    l_atom_pocket2 = parsePDB.loadCoordSectionPDB(path_file_pocket2) 
    
    l_atom_asa_pocket1 = parseNACCESS.fileASA(path_file_pocket1_asa)
    l_atom_asa_pocket2 = parseNACCESS.fileASA(path_file_pocket2_asa)    
    
    l_commun1 = []
    for dico_atom_pocket1 in l_atom_pocket1 : 
        for dico_atom_pocket2 in l_atom_pocket2 : 
            if dico_atom_pocket1["x"] == dico_atom_pocket2["x"] and dico_atom_pocket1["y"] == dico_atom_pocket2["y"] and dico_atom_pocket1["z"] == dico_atom_pocket2["z"] : 
                l_commun1.append (dico_atom_pocket1)
#                 l_commun2.append (dico_atom_pocket2)
                break
   
    score_g = 0.0
    score_inter = 0.0
    for dico_atom_asa in l_atom_asa_pocket1 : 
        score_g = score_g + dico_atom_asa["ABS"]
        for commun_atom in l_commun1 :
            if  dico_atom_asa["atomSeq"] == commun_atom["serial"] : 
                score_inter = score_inter + dico_atom_asa["ABS"]
                # break in case of same number residues
                break
    
    
    if debug :
        print len (l_commun1)
        print score_g
        print score_inter
        try : print score_inter / score_g
        except : pass
    
    
    if score_g == 0.0 : 
        return 0.0
    return score_inter / score_g
    


#MO ("/home/borrel/poche_lelieEstimator/SLE_pocket/1IVD_pocket11.pdb", "/home/borrel/poche_lelieEstimator/prox_sel/1IVD_pocket-NAG_472_NAG_473_A-NAG_484_NAG_485_A_atm.pdb","/home/borrel/poche_lelieEstimator/SLE_pocket/1IVD_pocket11_ACC.asa","/home/borrel/poche_lelieEstimator/prox_sel/1IVD_pocket-NAG_472_NAG_473_A-NAG_484_NAG_485_A_atm_ACC.asa")
     
#MO ("/home/borrel/poche_lelieEstimator/prox_sel/1IVD_pocket-NAG_472_NAG_473_A-NAG_484_NAG_485_A_atm.pdb", "/home/borrel/poche_lelieEstimator/SLE_pocket/1IVD_pocket11.pdb", "/home/borrel/poche_lelieEstimator/prox_sel/1IVD_pocket-NAG_472_NAG_473_A-NAG_484_NAG_485_A_atm_ACC.asa","/home/borrel/poche_lelieEstimator/SLE_pocket/1IVD_pocket11_ACC.asa")
  



   


from os import listdir, remove, system
from re import search

import runOtherProg
import getResidues
import overlapPocket


def retainPDB (listfile):
    listpdb = []
    for i in listfile:
        if i[-4:]==".pdb":
            listpdb.append(i)
    return(listpdb)


#~/Research/Projects/Topoisomerases/work/fev2017/pocketFiles/FullSet_pockets$ python ~/Research/Projects/src/srcOverlapScore/main_Overlap_arg.py pocketFiles/pockets_all/ pocketFiles/pockets_atp/ topo_GHKL_superpos_lig_ssATPLid_ssboucle80/ overlap/

def mainOverlap (p_dir_pocketEst1, p_dir_pocketEst2, p_dir_pro, p_result):
    

    
    # run NACCESS with complexe and pockets
    l_p_complexe = listdir(p_dir_pro)
    for file_pro in l_p_complexe : 
        if not search (".pdb", file_pro) : 
            continue
        # PDB ID
        PDB_ID = file_pro[0:4]
        print PDB_ID
        
        # open file overlapping
        filout_RO = open (p_result + file_pro[0:-4] + "_RO.txt", "w")
        filout_SO = open (p_result + file_pro[0:-4] + "_SO.txt", "w")
        filout_MO = open (p_result + file_pro[0:-4] + "_MO.txt", "w")
        
        
        p_pro = p_dir_pro + file_pro
        print p_pro
        
        # NACESS -> on complexe
        p_pro_naccess = runOtherProg.runNACESS(p_pro)
        
        # pocke est1
        l_p_pock1 = retainPDB (listdir(p_dir_pocketEst1))
        l_p_pock2 = retainPDB (listdir(p_dir_pocketEst2))
        
        # poke est2
        l_temp_file_pock1 = []
        for file_pock1 in l_p_pock1 :
            if search (PDB_ID, file_pock1) and search (".pdb", file_pock1): 
                p_files_pocket1_ACC  =  getResidues.getAccAtom(p_pro_naccess[0], p_dir_pocketEst1 + file_pock1)
                remove(p_files_pocket1_ACC[-1])
                l_temp_file_pock1.append (p_dir_pocketEst1 + file_pock1) 
        print l_temp_file_pock1    
        # pocke est2
        l_temp_file_pock2 = []
        for file_pock2 in l_p_pock2 : 
            if search (PDB_ID, file_pock2) and search (".pdb", file_pock2): 
                p_files_pocket2_ACC  =  getResidues.getAccAtom(p_pro_naccess[0], p_dir_pocketEst2 + file_pock2)
                remove(p_files_pocket2_ACC[-1])
                l_temp_file_pock2.append (p_dir_pocketEst2 + file_pock2) 
        
        
        
        # overlapping pocket
        #l_MO = []
        #l_RO = []
        #l_SO = []
        #flag_header_RO = 0
        
        # header file
        #if flag_header_RO == 0 :
        header = []
        for name_file_pock2 in l_temp_file_pock2 :
            header.append (name_file_pock2.split ("/")[-1][0:-4]) 
        filout_SO.write ("\t".join(header) + "\n")
        filout_RO.write ("\t".join(header) + "\n")    
        filout_MO.write ("\t".join(header) + "\n")
        
        
        for p_pocket_est1 in l_temp_file_pock1 : 
            l_MO = []
            l_SO = []
            l_RO = []
            for p_pocket_est2 in l_temp_file_pock2 :
                
                print "###### DEBUG ########"
                print PDB_ID
                print p_pocket_est1
                print p_pocket_est2
                print p_pocket_est1[0:-4] + "_ACC.asa"
                print p_pocket_est2[0:-4] + "_ACC.asa"        
                l_SO.append(str(overlapPocket.scoreOverlap(p_pocket_est1, p_pocket_est2)))
                l_MO.append(str(overlapPocket.MO(p_pocket_est1, p_pocket_est2, p_pocket_est1[0:-4] + "_ACC.asa", p_pocket_est2[0:-4] + "_ACC.asa")))
                l_RO.append (str(overlapPocket.MO(p_pocket_est2, p_pocket_est1, p_pocket_est2[0:-4] + "_ACC.asa", p_pocket_est1[0:-4] + "_ACC.asa")))
        
            filout_SO.write (p_pocket_est1.split ("/")[-1][0:-4] + "\t" + "\t".join(l_SO) + "\n")
            filout_MO.write (p_pocket_est1.split ("/")[-1][0:-4] + "\t" + "\t".join(l_MO) + "\n")
            print l_RO, "L_RO"
            print l_MO, "L_MO"
            print l_SO, "L_SO"
            filout_RO.write (p_pocket_est1.split ("/")[-1][0:-4] + "\t" + "\t".join(l_RO) + "\n")        
                
    
        filout_SO.close ()
        filout_MO.close ()
        filout_RO.close ()


def removeFileASA(path):
    cmdLR = "rm -rf "+path+"/*.asa " 
    system(cmdLR)
    cmdLR = "rm -rf "+path+"/*.rsa " 
    system(cmdLR)


# ----RUN----#
##############

#path_pocketSLE="/home/leslieregad/Research/Projects/EstimePocket/data/SLE_pockets/results_06-12-13/SL/"
#path_pocketFpocket="/home/leslieregad/Research/Projects/EstimePocket/data/fpocketProx_sel/"
#path_pocketProx="/home/leslieregad/Research/Projects/EstimePocket/data/prox_sel/"

#p_dir_pro="/home/leslieregad/Research/Projects/EstimePocket/data/structures/"
#path_res="/home/leslieregad/Research/Projects/EstimePocket/data/CompareEstime/OverlapScore/"



###overlap des poches Holo

#removeFileASA(p_dir_pro)
#mainOverlap (path_pocketSLE, path_pocketProx, p_dir_pro, path_res+"SLE_prox/")
#removeFileASA(p_dir_pro)
#mainOverlap (path_pocketProx, path_pocketFpocket, p_dir_pro, path_res+"Fpocket_prox/")
#removeFileASA(p_dir_pro)
#mainOverlap (path_pocketSLE, path_pocketFpocket, p_dir_pro, path_res+"SLE_Fpocket/")







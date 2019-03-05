#!/usr/bin/python


import sys


#pathsrc="/home/lregad/Research/Projects/src/srcOverlapScore/pocketgeneral/"

pathsrc="/home/hecini/Research/stage_HECINI/script/srcOverlapScore/pocketgeneral/"

sys.path.append(pathsrc)

from mainPocketOverlap_uniqSO_correspondingPocket import *

#import mainPocketOverlap_uniqSO_correspondingPocket

#----------------------------------#
#           personal path          #
#----------------------------------#

path_pocketSLE=sys.argv[1]   #directory qui contient les poches  estimated using the first approach
path_pocketProx=sys.argv[2]   #directory qui contient les poches estimated using the second approach
p_dir_pro=sys.argv[3] #directory containing the pdb files of proteins
path_res=sys.argv[4] #directory where the output files will be stored


#path_pocketSLE="/home/leslieregad/Research/Projects/CompareClassifPockLig/data/pocketsFiles/proximity/"
#path_pocketProx="/home/leslieregad/Research/Projects/CompareClassifPockLig/data/pocketsFiles/sle/"
#p_dir_pro="/home/leslieregad/Research/Projects/CompareClassifPockLig/data/proteinFilesOverLap/"
#path_res="/home/leslieregad/Research/Projects/CompareClassifPockLig/data/overlapScoreProx_SLE/"

#path_pocketSLE="pocketFiles/pockets_all/"
#path_pocketProx="pocketFiles/pockets_atp"
#p_dir_pro="pocketFiles/"
#path_res="overlap/"



#----------------------------------#
#               main               #
#----------------------------------#

removeFileASA(p_dir_pro)
mainOverlap (path_pocketSLE, path_pocketProx, p_dir_pro, path_res)




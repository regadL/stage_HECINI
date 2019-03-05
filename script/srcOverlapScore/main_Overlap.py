#!/usr/bin/python


import sys


pathsrc="/home/leslieregad/Research/Projects/src/srcOverlapScore/pocketgeneral/"
sys.path.append(pathsrc)
from mainPocketOverlap import *


#----------------------------------#
#           personal path          #
#----------------------------------#

path_pocketSLE="/home/leslieregad/Research/Projects/CompareClassifPockLig/data/pocketsFiles/proximity/"
path_pocketProx="/home/leslieregad/Research/Projects/CompareClassifPockLig/data/pocketsFiles/sle/"
p_dir_pro="/home/leslieregad/Research/Projects/CompareClassifPockLig/data/proteinFilesOverLap/"
path_res="/home/leslieregad/Research/Projects/CompareClassifPockLig/data/overlapScoreProx_SLE/"


#----------------------------------#
#               main               #
#----------------------------------#

removeFileASA(p_dir_pro)
mainOverlap (path_pocketSLE, path_pocketProx, p_dir_pro, path_res)




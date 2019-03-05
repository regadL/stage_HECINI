"""
BORREL Alexandre
8-10-12
main to run descriptors
"""
# personnal modules
import getResidues
import descriptor
import runOtherProg
import writeFiles
import loadDescriptors

# global modules
from os import system, path, makedirs, listdir
from re import search
from numpy.ma.core import indices




def main (path_file_pocket, path_complexe = "none", path_file_protomol = "none", name_ligand = "none", Renv = "none"):

    # path filout
    path_filout = path_file_pocket + ".desc"
    
    # temporaly directory
    path_dir_temp = path.dirname(path.abspath(path_filout)) + "/temp/"
    
    print path_dir_temp
    try :makedirs( path_dir_temp, mode=0777 )
    except :pass



    # run NACCESS
    if path_complexe != "none" : 
        path_files_naccess = runOtherProg.runNACESS(path_complexe)
    else :
        path_files_naccess = runOtherProg.runNACESS(path_file_pocket)
    
    
    path_file_pocket_ACC  =  getResidues.getAccAtom(path_files_naccess[0], path_file_pocket)
#    
#    
    descriptor.runDescriptor(path_file_pocket, path_file_pocket, path_file_pocket_ACC[1], path_filout, path_complexe, path_file_protomol, name_ligand)

    system ("mv " + path_file_pocket_ACC[0] + " " + path_dir_temp + path.basename(path_file_pocket_ACC[0]) ) 
    system ("mv " + path_file_pocket_ACC[1] + " " + path_dir_temp + path.basename(path_file_pocket_ACC[1]) )
    system ("mv " + path_files_naccess[0] + " " + path_dir_temp + path.basename(path_files_naccess[0]) ) 
    system ("mv " + path_files_naccess[1] + " " + path_dir_temp + path.basename(path_files_naccess[1]) )
    
    
    # run druggability model
    if Renv != "none" : 
        p_descriptor = path_file_pocket + ".Rdesc"
        dico_desc = {} 
        loadDescriptors.loadDescriptorPocket (path_desc_file = path_file_pocket + ".desc", dico_load = dico_desc)
        writeFiles.onePocketDescriptor (dico_desc[""], p_descriptor)
        
        runOtherProg.predictLDA(Renv, p_descriptor, path_file_pocket + ".predict")



#path_dir = "/home/borrel/Dropbox/TeamDruggability/dataLeslie/"

# main ("/home/borrel/Desktop/pocket0_atm.pdb", path_complexe = "none", path_file_protomol = "none", name_ligand = "none")

# main ("/home/borrel/druggabilityProject/modelApplication/pocket4/1R58_A41_pocket.pdb", path_complexe = "/home/borrel/druggabilityProject/modelApplication/pocket4/protein/1R58.pdb", Renv = "/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/WithLigand/with_ligand.Rdata")
# main ("/home/borrel/druggabilityProject/modelApplication/pocket4/1R58_AO5_pocket.pdb", path_complexe = "/home/borrel/druggabilityProject/modelApplication/pocket4/protein/1R58.pdb", Renv = "/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/WithLigand/with_ligand.Rdata")
# main ("/home/borrel/druggabilityProject/modelApplication/pocket4/1R58_TN4_pocket.pdb", path_complexe = "/home/borrel/druggabilityProject/modelApplication/pocket4/protein/1R58.pdb", Renv = "/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/WithLigand/with_ligand.Rdata")
# main ("/home/borrel/druggabilityProject/modelApplication/pocket4/1YW7_superpos1R58_A41_pocket.pdb", path_complexe = "/home/borrel/druggabilityProject/modelApplication/pocket4/protein/1YW7_superpos1R58.pdb", Renv = "/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/WithLigand/with_ligand.Rdata")
# main ("/home/borrel/druggabilityProject/modelApplication/pocket4/3FMR_superpos1R58_TN4_pocket.pdb", path_complexe = "/home/borrel/druggabilityProject/modelApplication/pocket4/protein/3FMR_superpos1R58.pdb", Renv = "/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/WithLigand/with_ligand.Rdata")

# list_files = listdir (path_dir + "Pocket/")
# print list_files
# 
# print list_files.index ("AnaB-AHBD_B_FAD_Allpockets_ATOM.pdb")
# 
# print len (list_files)
# 
# for file_dir in list_files[40:43] : 
#     if search ("ATOM.pdb$", file_dir) : 
#         print file_dir
#         main (path_dir + "Pocket/"+ file_dir, path_complexe = path_dir + "Prot/" + file_dir.split ("_")[0] + ".pdb")




path_dir = "/home/borrel/Bureau/PourAlex/"
pdb_dir = path_dir + "pdbFile/"
pocket_dir = path_dir + "pocketFile/"

l_file_pocket = listdir(pocket_dir)
l_file_pr = listdir (pdb_dir)



# for file_pocket in l_file_pocket : 
#     pdb = file_pocket.split ("_")[0]
#     try : main (pocket_dir + file_pocket, path_complexe = pdb_dir + pdb + "_A.pdb", Renv = "/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/WithoutLigand/without_ligand.Rdata")
#     except : main (pocket_dir + file_pocket, path_complexe = pdb_dir + pdb + "_B.pdb", Renv = "/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/WithoutLigand/without_ligand.Rdata")

filout = open (path_dir + "result", "w")

for file_pocket in l_file_pocket : 
    if search(".predict", file_pocket) : 
        filin = open (pocket_dir + file_pocket, "r")
        l_lines = filin.readlines ()
        filin.close ()
        
        proba = l_lines[-1].split ()[2]
        print proba
        filout.write (file_pocket + " " + str (proba) + "\n")

filout.close ()

    






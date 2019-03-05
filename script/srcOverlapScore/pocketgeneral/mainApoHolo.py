"""
BORREL Alexandre
09-2012
"""

import globalFonction
import runOtherProg
import writeFiles
import pathDirectory
import tool
import analysis
import superposeStructure

from os import listdir, path, remove
from re import search

def main (path_file_model, name_dataset="ApoForm"):
    
    
    dico_dataset = globalFonction.calculDatasetDictionary(name_dataset,1)
    
#    print len(dico_dataset.keys())

    path_file_correspondace = writeFiles.corespondanceApoHolo (dico_dataset, pathDirectory.result(name_dataset) + "correspondencePDB")
    
    # divise dataset    
    
    dico_dataset_type_apo = tool.selectOnlyTypeStructure(dico_dataset, "apo structure")
    dico_dataset_type_holo = tool.selectOnlyTypeStructure(dico_dataset, "holo structure") 
    
    print len (dico_dataset_type_apo.keys ())
    print dico_dataset_type_holo[dico_dataset_type_holo.keys ()[1]]
        
    
    # accessibility solvent
    runOtherProg.globalNACCESS(dico_dataset_type_apo, name_dataset)
    runOtherProg.globalNACCESS(dico_dataset_type_holo, name_dataset)
    
    # superimpose -> retrieve matrix transloc
    superposeStructure.superposeApoHolo(dico_dataset, name_dataset)
    
    ######################
    # Pocket Estimation  #
    ######################    
    
    # pocket estimation holo
    globalFonction.pocketEstimation(name_dataset, dico_dataset_type_holo, "Fpocket", runFpocket = 1)
    
    # pocket estimation apo
    globalFonction.pocketEstimationApoForm(name_dataset, dico_dataset_type_apo, dico_dataset_type_holo, runFpocket = 1)
    
    
    ##############
    # Descriptor #
    ##############
    
    # descriptor
    dico_descriptors_type_apo = globalFonction.retrieveGlobalDescriptors ("Fpocket", "none", dico_dataset_type_apo, name_dataset, write_file = 0, calcul_descriptor = 1, option_separate = 1)
    dico_descriptors_type_holo = globalFonction.retrieveGlobalDescriptors ("Fpocket", "none", dico_dataset_type_holo, name_dataset, write_file = 0, calcul_descriptor = 1, option_separate = 1)

    print dico_descriptors_type_holo

    #####################
    # write data global #
    #####################
    path_file_descriptor_apo = writeFiles.globalDescriptors(dico_descriptors_type_apo, pathDirectory.result(name_dataset) + "apo_all_pocket.data")
    path_file_descriptor_holo = writeFiles.globalDescriptors(dico_descriptors_type_holo, pathDirectory.result(name_dataset) + "holo_all_pocket.data")
  
    path_file_RMSD = pathDirectory.searchRMSDFile (pathDirectory.descriptor(name_dataset))
    
    path_dir_result = pathDirectory.result(name_dataset)
    # color for ACP

    writeFiles.colorACPFile (dico_descriptors_type_apo)
       
    ###############
    #     ACP     #
    ###############
###    # apo protein
    analysis.specificACP("global", dico_descriptors_type_apo, path_dir_result + "apo_globalData", mainACP = "All_descriptors_Apo")
    analysis.specificACP(["hydrophobic_kyte","p_Ooh_atom", "p_aromatic_residues"], dico_descriptors_type_apo, path_dir_result + "apo_descModel", mainACP = "Specific_descriptors")
   
###    # holo protein
    # analysis.specificACP("global", dico_descriptors_type_holo, path_dir_result + "holo_globalData", mainACP = "All_descriptors_holo")
    # analysis.specificACP(["hydrophobic_kyte","p_Ooh_atom", "p_aromatic_residues"], dico_descriptors_type_holo, path_dir_result + "holo_descModel", mainACP = "Specific_descriptors")
   
###    # PCA apo and holo same plot
    # runOtherProg.ACPDataset (path_file_descriptor_apo, path_file_descriptor_holo, path_dir_result + "PCA_apo_holo")
    analysis.ACPTwoDatasetDescriptor(dico_descriptors_type_apo, dico_descriptors_type_holo, ["hydrophobic_kyte","p_Ooh_atom", "p_aromatic_residues"], path_dir_result + "desc_model", correspondance_file=path_file_correspondace)
    analysis.ACPTwoDatasetDescriptor(dico_descriptors_type_apo, dico_descriptors_type_holo, "radi", path_dir_result+ "radi", correspondance_file=path_file_correspondace)
    analysis.ACPTwoDatasetDescriptor(dico_descriptors_type_apo, dico_descriptors_type_holo,["RADIUS_HULL", "DIAMETER_HULL", "SURFACE_HULL", "VOLUME_HULL", "SMALLEST_SIZE", "INERTIA_3", "INERTIA_1", "FACE", "PCI", "PSI", "RADIUS_CYLINDER", "X._ATOM_CONVEXE", "CONVEX.SHAPE_COEFFICIENT", "INERTIA_2", "C_RESIDUES", "C_ATOM"], path_dir_result+ "geoOnly", correspondance_file=path_file_correspondace)
#      
    
    
    ################
    # RMSD pocket  #
    ################
    
    analysis.RMSDPockets (dico_dataset, pathDirectory.result(name_dataset + "/RMSD"), name_dataset)
    
   
    ########################
    # Histogram descriptor #
    ########################
    # analysis.histogramFonctionRMSD (dico_descriptors_type_apo, dico_descriptors_type_holo, path_file_RMSD, path_dir_result, ["hydrophobic_kyte"])
    # analysis.histogramFonctionRMSD (dico_descriptors_type_apo, dico_descriptors_type_holo, path_file_RMSD, path_dir_result, ["p_Ooh_atom"])
    # analysis.histogramFonctionRMSD (dico_descriptors_type_apo, dico_descriptors_type_holo, path_file_RMSD, path_dir_result, ["p_aromatic_residues"])
    
   
    ###################
    # apply LDA model #
    ###################
#    # global
#    # apo
     
    globalFonction.applyModel(path_file_model, path_file_descriptor_apo, path_dir_result + "predictApo.result") # laisse seulement une commande dans la fonction pour apres car surement ACP ou autre
#    
#    # holo
    globalFonction.applyModel(path_file_model, path_file_descriptor_holo, path_dir_result + "predictHolo.result")
   #by type holo structure
    
 
 
 
    
    
def multiModelTest (path_dir_model, name_file_result = "", name_dataset="ApoForm"): 
    
    f = 0
    dico_dataset = globalFonction.calculDatasetDictionary(name_dataset,0)
    if name_dataset == "ApoForm138" : 
        name_dataset = "ApoForm"
        f = 1
    
    path_file_correspondace = writeFiles.corespondanceApoHolo (dico_dataset, pathDirectory.result(name_dataset) + "correspondencePDB")
    
    # divise dataset    
    dico_dataset_type_apo = tool.selectOnlyTypeStructure(dico_dataset, "apo structure")
    dico_dataset_type_holo = tool.selectOnlyTypeStructure(dico_dataset, "holo structure") 
    
    # dictionnary with descriptors   
    dico_descriptors_type_apo = globalFonction.retrieveGlobalDescriptors ("Fpocket", "none", dico_dataset_type_apo, name_dataset, write_file = 0, calcul_descriptor = 0, option_separate = 1)
    dico_descriptors_type_holo = globalFonction.retrieveGlobalDescriptors ("Fpocket", "none", dico_dataset_type_holo, name_dataset, write_file = 0, calcul_descriptor = 0, option_separate = 1)

    if f == 1 : 
        name_dataset = "ApoForm138"

    path_dir_result = pathDirectory.result(name_dataset)
    #####################
    # write data global #
    #####################
    path_file_descriptor_apo = writeFiles.globalDescriptors(dico_descriptors_type_apo, pathDirectory.result(name_dataset) + "apo_all_pocket.data")
    path_file_descriptor_holo = writeFiles.globalDescriptors(dico_descriptors_type_holo, pathDirectory.result(name_dataset) + "holo_all_pocket.data")
  

    l_file_model = listdir(path_dir_model)
    p_file_result = path_dir_result + "best" + str (name_file_result) + ".result"
    # check exist ?
    if path.exists(p_file_result) : 
        remove(p_file_result)

    for p_file_model in l_file_model : 
        if search("Rdata", p_file_model) : 
            runOtherProg.predictLDA (path_dir_model + p_file_model, path_file_descriptor_apo, p_file_result, plot = 0)


# Dataset Apo Schmitke global -> refaire la fonction de codage  #
#################################################################

# main ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/AutoSelected/autoselected.Rdata", name_dataset="ApoForm")
# main ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/AutoSelected/autoselected.Rdata", name_dataset="ApoForm138")


# Dataset Apo Schmitke clean -> only one apo by holo
# main ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/AutoSelected/autoselected.Rdata", name_dataset="ApoFormClean")
#####main ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/hydroAro/hydroAro.Rdata", name_dataset="ApoFormClean")

# apo huang
# main ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/AutoSelected/autoselected.Rdata", name_dataset = "ApoHuang")


##############################################
#  test best models proposed with selection  #
##############################################
# -POE- #
# multiModelTest ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/selectedDesc/BestModels25/", name_dataset="ApoForm138", name_file_result = "POE3_25")
# multiModelTest ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/selectedDesc4/BestModels25/", name_dataset="ApoForm138", name_file_result = "POE4_25")

# multiModelTest ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/selectedDesc/BestModels25/", name_dataset="ApoFormClean", name_file_result = "POE3_25")
# multiModelTest ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/selectedDesc4/BestModels25/", name_dataset="ApoFormClean", name_file_result = "POE4_25")

# -PLE- #
# multiModelTest ("/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/selectedDesc/BestModels25/", name_dataset="ApoForm138", name_file_result = "PLE_25")
# multiModelTest ("/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/selectedDesc/BestModels25/", name_dataset="ApoFormClean", name_file_result = "PLE_25")

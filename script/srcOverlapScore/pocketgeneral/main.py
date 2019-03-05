# personal module
import globalFonction
import analysis
import runOtherProg
import tool
import loadDescriptors
import comparison
import rateQuality
import pathDirectory
import ligandDrugLike
import zipCompress



#################################
#     Fpocket parameters        #
#################################

#analysis.fpocketParameters(dico_K.keys ()[0:50], run = 0, step = 0.1)

#### revoir la selection de variable pour plus performant


def main (name_dataset, protomol_type, pocket_retrive_type) : 
    
    #################
    # Download data #
    #################
    # download file and generate dataset
    dico_dataset = globalFonction.calculDatasetDictionary(name_dataset,0)
#     print dico_dataset
#    print dico_dataset
#    i = 0
#    print len (dico_dataset.keys ())
#     for PDB in dico_dataset.keys (): 
#         print PDB, dico_dataset[PDB]
#        if dico_dataset[PDB]["Type structure"] == "apo structure" : 
#            print PDB, dico_dataset[PDB]
#            i = i +1
#    print i
#        
#    print list_PDB.index("1Z8A") # for test
    #######################
    #  Analysis Data Set  #
    #######################
#    globalFonction.analysisDataset(dico_dataset, name_dataset)
    
    ######################
    # Pocket Estimation  #
    ######################
#     globalFonction.pocketEstimation(name_dataset, dico_dataset, pocket_retrive_type, runFpocket = 1)

    ############################
    # Run accesibility solvant #
    ############################
#    runOtherProg.globalNACCESS(dico_dataset, name_dataset)

    #################################
    # Run Surflex generate protomol #
    #################################
##    globalFonction.generationProtomol (dico_dataset, pocket_retrive_type, name_dataset)


    ###############################
    #  Analysis pocket estimation #
    ###############################    
#     globalFonction.analysisPocketEstimation(dico_dataset, name_dataset, pocket_retrive_type)
    
    ###############################
    #      RUN Descriptors        #
    ###############################
    dico_descriptors = globalFonction.retrieveGlobalDescriptors (pocket_retrive_type, protomol_type, dico_dataset, name_dataset, write_file = 0, calcul_descriptor = 0, option_separate = 1, compo_aa = 0)
    
#     print dico_descriptors.keys ()
#     print dico_descriptors["data"]
#     print dico_descriptors["Druggable"]["PDB"].index ("1D09")
#     print dico_descriptors["Druggable"]["data"]
    
#     print dico_descriptors["Druggable"]["data"][62], dico_descriptors["Druggable"]["PDB"][62]
    #del rugosity (incomplet)
#     try :
#         del dico_descriptors["Druggable"]["area"]["rugosity"]
#         del dico_descriptors["No-Druggable"]["area"]["rugosity"]
#     except : 
#         pass
##    
#    ########################
#    #      Analysis        #
#    ########################
    # run prediction perhaps make new function to separate visualization / prediction / correlation
    globalFonction.analysisGlobalDescriptor (dico_descriptors, name_dataset, pocket_retrive_type, ttest = 0, ACP = 0, LDAGlobal = 1, CART = 0, correlation = 0, selection_variable = 1, SVM = 0, histogram = 0, radomForest = 0, regLog = 0)
    

    #####################################
    #   Specific descriptors analysis   #
    #####################################
#    analysis.correlationVolumeVSVolumeFpocket(dico_descriptors, dico_dataset.keys (), name_dataset, pocket_retrive_type)
    

def model (name_dataset_train, name_dataset_test, pocket_retrieve_type, protomol_type, pvalue, cor_value):
    """
    A voir pour changer et avoir un truc plus propre
    """
    globalFonction.validation(name_dataset_train, name_dataset_test, pocket_retrieve_type, protomol_type, pvalue, cor_value)
    
   

#######################

data_set = "krasowski"
#main (data_set, "surflexe", "proximity")
#main (data_set, "surflexe", "Fpocket")
#
main (data_set, "none", "proximity")
main (data_set, "none", "Fpocket")
# main (data_set, "none", "cavitator")

#Case wih dogsite -> unzio the pocket folder
# zipCompress.unCompressFolder("/home/borrel/druggabilityProject/DOGSITE/")
# main (data_set, "none", "DogSite")

data_set = "DD"
# main (data_set, "none", "Fpocket")
# main (data_set, "none", "proximity")

data_set = "Perola"

#main (data_set, "none", "Fpocket")
#main (data_set, "surflexe", "proximity")


######################
# DRUGGABILITY MODEL #
######################

data_train = "krasowski"
data_test = "DD"
pocket_type = "proximity"
protomol = "none"

#model (data_train, data_test, pocket_type, protomol, 0.10, 1)
#model (data_train, data_train, pocket_type, protomol, 0.10, 1)

pocket_type = "Fpocket"
#model (data_train, data_test, pocket_type, protomol, 0.10, 1)
#model (data_train, data_train, pocket_type, protomol, 0.10, 1)

data_train = "krasowski_train"
data_test = "krasowski_test"

#model (data_train, data_test, pocket_type, protomol, 0.10, 1)


# Attention -> ajouter un filtre pour suprimer les descripteurs qui ont des NA (juste la rugosite quand le plan ne peut pas etre calculer par manque de points) a revoir

###############cross analysis###############
############################################

# ACP with 2 dataset -> 16-10 change fonction
# descriptor selected
#globalFonction.ACPDataset1BetwwenDataset2 ("krasowski", "krasowski", "proximity", "none", "Fpocket", "none", list_descriptor = ["p_Ooh_atom", "hydrophobic_kyte", "RADIUS_CYLINDER", "PCI", "p_aromatic_residues"])
#globalFonction.ACPDataset1BetwwenDataset2 ("krasowski", "krasowski", "proximity", "none", "Fpocket", "none", list_descriptor = ["RADIUS_HULL", "DIAMETER_HULL", "SURFACE_HULL", "VOLUME_HULL", "SMALLEST_SIZE", "INERTIA_3", "INERTIA_1", "FACE", "HEIGHT_CYLINDER", "PCI", "PSI", "RADIUS_CYLINDER", "%_ATOM_CONVEXE", "CONVEX-SHAPE_COEFFICIENT", "INERTIA_2", "c_residues", "c_atom", "Real_volume"])
#globalFonction.ACPDataset1BetwwenDataset2 ("krasowski", "krasowski", "proximity", "none", "Fpocket", "none", list_descriptor = ["RADIUS_HULL", "DIAMETER_HULL", "SURFACE_HULL", "VOLUME_HULL", "SMALLEST_SIZE", "INERTIA_3", "INERTIA_1", "FACE", "HEIGHT_CYLINDER", "PCI", "PSI", "RADIUS_CYLINDER", "%_ATOM_CONVEXE", "CONVEX-SHAPE_COEFFICIENT", "INERTIA_2", "c_residues", "c_atom"])


# All descriptors
# globalFonction.ACPDataset1BetwwenDataset2 ("krasowski", "krasowski", "proximity", "none", "Fpocket", "none")
# globalFonction.MDSDataset1BetwwenDataset2 ("krasowski", "krasowski", "proximity", "none", "Fpocket", "none")

# globalFonction.ACPDataset1BetwwenDataset2 ("krasowski", "krasowski", "proximity", "none", "Fpocket", "none", list_descriptor=["hydrophobic_kyte", "p_hydrophobic_residues" ,"p_hydrophobic_atom", "p_hyd_atom","hydrophobicity_pocket_pocket", "p_aromatic_residues", 
#                               "p_Car_atom", "p_polar_residues" ,"p_aliphatic_residues","p_Nlys_atom", "p_Ntrp_atom", "p_S_atom", "p_Otyr_atom", "p_Ooh_atom",
#                                "p_O_atom", "p_N_atom", "p_ND1_atom", "p_NE2_atom", "polarity_pocket_pocket", "p_charged_residues", "p_positive_residues", "p_negative_residues",
#                                "p_Ocoo_atom", "p_Cgln_atom",  "p_Ccoo_atom", "p_Carg_atom", "charge",  "p_pro_residues", "p_tiny_residues",  "p_main_chain_atom", 
#                                "p_side_chain_atom", "P_C_atom", "p_nitrogen_atom", "p_sulfur_atom", "p_oxygen_atom", "p_carbone_atom"])




# globalFonction.MDSDataset1BetwwenDataset2 ("krasowski", "krasowski", "proximity", "none", "Fpocket", "none")


#globalFonction.ACPDataset1BetwwenDataset2 ("krasowski", "Perola", "proximity", "surflexe") !!!!!
#globalFonction.ACPDataset1BetwwenDataset2 ("krasowski", "Perola", "proximity", "surflexe") !!!!!!
#globalFonction.ACPDataset1BetwwenDataset2 ("krasowski", "DD", "proximity", "none") !!!!!!!





#################
## Correlation ##
#################

# correlation with 2 types estimation
#globalFonction.correlationDescriptorsByTypePocket ("krasowski", "Fpocket", "surflexe",  "proximity", "surflexe" )
#
# globalFonction.correlationDescriptorsByTypePocket ("krasowski", "Fpocket", "none",  "proximity", "none" )



### Proximity ####
##################
#globalFonction.correlationDescriptor("krasowski", ["hydrophobic_kyte", "aromatic_residues"],"proximity", "none")
#globalFonction.correlationDescriptor("krasowski", ["hydrophobic_kyte", "PCI"],"proximity", "none")
#globalFonction.correlationDescriptor("krasowski", ["hydrophobic_kyte", "RADIUS_CYLINDER"],"proximity", "none")
##
#globalFonction.correlationDescriptor("krasowski", ["RADIUS_CYLINDER", "PCI"],"proximity", "none")
#globalFonction.correlationDescriptor("krasowski", ["aromatic_residues", "RADIUS_CYLINDER"],"proximity", "none")

### Fpocket ###
###############
#globalFonction.correlationDescriptor("krasowski", ["hydrophobic_kyte", "aromatic_residues"],"Fpocket", "none")
#globalFonction.correlationDescriptor("krasowski", ["hydrophobic_kyte", "alcool"],"Fpocket", "none")


##################
# FPOCKET SCORE #
#################

#rateQuality.qualityFpocket("Fpocket", "krasowski", pathDirectory.result("Score_Fpocket") + "quality_Fpocket_krasowski")
# rateQuality.qualityFpocket("Fpocket", "DD", pathDirectory.result("comparison_with_Fpocket") + "quality_Fpocket_DD")

# rateQuality.qualityFpocket("Fpocket", "ApoForm138", pathDirectory.result("comparison_with_Fpocket") + "quality_Fpocket_apo138")


################
# comparison   #
################

# comparison.comparisonWithLigandWithoutLigand("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/WithoutLigand","/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/WithLigand")
# comparison.FpocketWithoutLigand ("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/WithoutLigand", "krasowski", "Fpocket", "none")


####################
#  LIPINSKI RULES  #
####################
#dico_dataset = globalFonction.calculDatasetDictionary("krasowski",0)
#path_file_list_ligand = ligandDrugLike.retrieveSmileDrugLike (dico_dataset, pathDirectory.dataSet("krasowski") + "ligandAnalysis") [1] # list files
# After weeb service http://crdd.osdd.net:8081/webcdk/desc.jsp
#path_file_list_ligand = "/home/borrel/druggabilityProject/dataSet/krasowski/ligandAnalysis.list"
#
#path_ligand_druglike = ligandDrugLike.parseLigandDescriptor (path_file_list_ligand, "/home/borrel/druggabilityProject/dataSet/krasowski/ligand_descripteur.txt", "/home/borrel/druggabilityProject/dataSet/krasowski/ligand_lipinski.txt")


#####################
#  OVERLAAP POCKET  #
#####################

# globalFonction.overlapTwoPockets("krasowski", "Fpocket", "proximity")



#########################
#  Accesibility Pocket  #
#########################

# globalFonction.retrieveAccessibilityPocket("krasowski", "proximity")


####################
#   APPLY MODEL    #
####################

# globalFonction.applyModel("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/best2/without_ligand.Rdata", "/home/borrel/druggabilityProject/result/DD/Fpocket/applyModel/global.data", "/home/borrel/druggabilityProject/result/DD/Fpocket/applyModel/best2predictFpocket.result", family_barplot=0, name_dataset="DD")
# globalFonction.applyModel("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/best3/without_ligand.Rdata", "/home/borrel/druggabilityProject/result/DD/Fpocket/applyModel/global.data", "/home/borrel/druggabilityProject/result/DD/Fpocket/applyModel/best3predictFpocket.result", family_barplot=0, name_dataset="DD")

# globalFonction.applyModel("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/WithoutLigand/without_ligand.Rdata", "/home/borrel/druggabilityProject/result/DD/Fpocket/applyModel/without_ligand", "/home/borrel/druggabilityProject/result/DD/Fpocket/applyModel/predictFpocket.result", family_barplot=1, name_dataset="DD")
# globalFonction.applyModel("/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/WithLigand/with_ligand.Rdata", "/home/borrel/druggabilityProject/result/DD/proximity/applyModel/with_ligand", "/home/borrel/druggabilityProject/result/DD/proximity/applyModel/predictProx.result")






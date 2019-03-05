"""
BORREL Alexandre
04-2012
"""
# personal modules
import analysis
import pathDirectory
import downloadFile
import tool
import checkPocket
import dataSet
import writeFiles
import preparePDB
import identityCalcul
import runOtherProg
import getResidues
import descriptor
import loadDescriptors
import parameters
import superposeStructure
import overlapPocket

# general module
import os, re
import descriptorEnergy


def analysisGlobalDescriptor (dictionary_descriptor, dataset_name, pocket_type, ttest = 1, ACP = 1, LDAGlobal = 1, CART = 1, SVM = 1,  correlation = 1, histogram = 1, selection_variable = 0, radomForest = 1, regLog = 1):
    
    path_file_color = writeFiles.colorACPFile (dictionary_descriptor) # color for ACP
    
    list_with_ligand = ["hydrophobic_kyte", "PCI", "RADIUS_CYLINDER", "p_aromatic_residues"]
    list_without_ligand = ["hydrophobic_kyte","p_Ooh_atom", "p_aromatic_residues"]
    
    
    
    if ttest : 
        path_result = pathDirectory.result (dataset_name + "/" + pocket_type + "/Ttest")
        if len (dictionary_descriptor["Druggable"].keys()) > 1 and len (dictionary_descriptor["No-Druggable"].keys()) > 1 :  
            dico_ttest = analysis.ttest(dictionary_descriptor["Druggable"], dictionary_descriptor["No-Druggable"], path_result + "ttest_" + dataset_name )
            list_result_ttest_signif = analysis.retrieveDescriptorSignifTtest (dico_ttest, pvalue=0.10)
            list_descriptor_ttest_signif = list_result_ttest_signif.keys () # because confidence interval is calculate
        print list_descriptor_ttest_signif
    
    
    if ACP : 
        path_result = pathDirectory.result (dataset_name + "/" + pocket_type + "/ACP")
        analysis.specificACP("global", dictionary_descriptor, path_result + "des_global", mainACP = "Descriptor_global")
        
#        if 'list_descriptor_ttest_signif' in locals () :
#            analysis.specificACP(list_descriptor_ttest_signif, dictionary_descriptor, path_result + "signifttest", mainACP = "T-test_signif")
#            
#        analysis.specificACP("volume", dictionary_descriptor, path_result + "des_volume", mainACP = "Descriptor_volume")
#        analysis.specificACP("atomic", dictionary_descriptor, path_result + "des_atomic", mainACP = "Descriptor_atomic")
#         analysis.specificACP("fpocket", dictionary_descriptor, path_result + "des_fpocket", mainACP = "Descriptor_Fpocket")  
#        analysis.specificACP(["RADIUS_HULL", "DIAMETER_HULL", "SURFACE_HULL", "VOLUME_HULL", "SMALLEST_SIZE", "INERTIA_3", "INERTIA_1", "FACE", "HEIGHT_CYLINDER", "PCI", "PSI", "RADIUS_CYLINDER", "%_ATOM_CONVEXE", "CONVEX-SHAPE_COEFFICIENT", "INERTIA_2", "c_residues", "c_atom", "Real_volume"], dictionary_descriptor, path_result + "des_geo", mainACP = "Descriptor_Fpocket") 
#         analysis.specificACP(["C_residues", "C_atom"], dictionary_descriptor, path_result + "des_count", mainACP = "Descriptor_count") 
        analysis.specificACP(["VOLUME_HULL", "SMALLEST_SIZE", "RADIUS_HULL", "DIAMETER_HULL" ,"SURFACE_HULL", "RADIUS_CYLINDER" , "C_ATOM", "C_RESIDUES", 
                              "PSI", "INERTIA_2", "INERTIA_3", "INERTIA_1", "CONVEX.SHAPE_COEFFICIENT", "PCI", "FACE", "X._ATOM_CONVEXE"],
                             dictionary_descriptor, path_result + "des_vol_form", mainACP = "Descriptor_vol_form")
        
        analysis.specificACP(["hydrophobic_kyte", "p_hydrophobic_residues" ,"p_hydrophobic_atom", "p_hyd_atom","hydrophobicity_pocket_pocket", "p_aromatic_residues", 
                              "p_Car_atom", "p_polar_residues" ,"p_aliphatic_residues","p_Nlys_atom", "p_Ntrp_atom", "p_S_atom", "p_Otyr_atom", "p_Ooh_atom",
                               "p_O_atom", "p_N_atom", "p_ND1_atom", "p_NE2_atom", "polarity_pocket_pocket", "p_charged_residues", "p_positive_residues", "p_negative_residues",
                               "p_Ocoo_atom", "p_Cgln_atom",  "p_Ccoo_atom", "p_Carg_atom", "charge",  "p_pro_residues", "p_tiny_residues",  "p_main_chain_atom", 
                               "p_side_chain_atom", "P_C_atom", "p_nitrogen_atom", "p_sulfur_atom", "p_oxygen_atom", "p_carbone_atom"], dictionary_descriptor, path_result + "des_physico", mainACP = "Descriptor_des_physico")




        if pocket_type == "Fpocket" : 
            analysis.specificACP (list_without_ligand, dictionary_descriptor,  path_result +"without_ligand", "descriptor_selected") 
        else : 
            analysis.specificACP (list_with_ligand,dictionary_descriptor, path_result +"with_ligand", "descriptor_selected") 
     
     
     
        
    if LDAGlobal :
        
        # path directory
        path_begin = pathDirectory.result (dataset_name + "/" + pocket_type + "/LDA")
        # global data for Rscript
        path_global_descriptor = writeFiles.globalDescriptors(dictionary_descriptor, path_begin + "global.data")
        
        if selection_variable == 1 : 
            path_dir_parameter = pathDirectory.generatePath(path_begin + "selectedDesc/")
            print path_dir_parameter
            
            list_path_files = writeFiles.specificDescriptorbyData(dictionary_descriptor, "global", path_dir_parameter + "global_desc")
            list_descriptor_selected = parameters.retrieveDescriptorForLDAModel (dictionary_descriptor, list_path_files[0], list_path_files[1], path_dir_parameter, nb_model = 10)   
        
#         analysis.LDAGlobal (dictionary_descriptor, "global", pathDirectory.generatePath(path_begin + "DescGlobal/"), "des_global", "Descriptor_global", test_train_data = 1, name_dataset = dataset_name, path_file_global_descriptor = path_global_descriptor)
#        analysis.LDAGlobal (dictionary_descriptor, "radi", pathDirectory.generatePath(path_begin + "DescRadi/"), "des_radi", "Descriptor_radi",test_train_data = 1, name_dataset = dataset_name, path_file_global_descriptor = path_global_descriptor)
#       # descriptor selected
#         if "list_descriptor_selected" in locals () :
#             analysis.LDAGlobal (dictionary_descriptor, list_descriptor_selected, pathDirectory.generatePath(path_begin + "AutoSelected/"), "autoselected", "autoselected", path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
#        
#        # Descriptors imposed
#         if pocket_type == "Fpocket" : 
#             analysis.LDAGlobal (dictionary_descriptor, list_without_ligand, pathDirectory.generatePath(path_begin + "WithoutLigand/"), "without_ligand", "without_ligand",  path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
#         else : 
#             analysis.LDAGlobal (dictionary_descriptor, list_with_ligand, pathDirectory.generatePath(path_begin + "WithLigand/"), "with_ligand", "with_ligand",  path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
       
#         analysis.LDAGlobal (dictionary_descriptor, ["hydrophobic_kyte","p_Ooh_atom","p_main_chain_atom"], pathDirectory.generatePath(path_begin + "best2/"), "without_ligand", "without_ligand",  path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
#         analysis.LDAGlobal (dictionary_descriptor, ["hydrophobic_kyte", "p_Carg_atom", "p_aromatic_residues"], pathDirectory.generatePath(path_begin + "best3/"), "without_ligand", "without_ligand",  path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)

       
        
        # test possibility with less descriptors
#        analysis.LDAGlobal (dictionary_descriptor, ["hydrophobic_kyte"], pathDirectory.generatePath(path_begin + "hydroOnly/"), "hydroOnly", "hydroOnly",  path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
#         analysis.LDAGlobal (dictionary_descriptor, ["hydrophobic_kyte", "p_aromatic_residues"], pathDirectory.generatePath(path_begin + "hydroAro/"), "hydroAro", "hydroAro",  path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
#         analysis.LDAGlobal (dictionary_descriptor, ["hydrophobic_kyte", "RADIUS_CYLINDER"], pathDirectory.generatePath(path_begin + "hydroRadius/"), "hydroRadius", "hydroRadius",  path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
#         analysis.LDAGlobal (dictionary_descriptor, ["hydrophobic_kyte", "PCI", "RADIUS_CYLINDER"], pathDirectory.generatePath(path_begin + "hydroRadiusPCI/"), "hydroRadiusPCI", "hydroRadiusPCI",  path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)


        
    if CART : 
        path_begin = dataset_name + "/" + pocket_type + "/CART/"
        analysis.cart(dictionary_descriptor, "global", pathDirectory.result (path_begin + "des_global") + "des_global")
        
#      
        if "list_descriptor_ttest_signif" in locals () :
            analysis.cart (dictionary_descriptor, list_descriptor_ttest_signif, pathDirectory.result (path_begin + "des_Ttest") + "des_Ttest")
            
        if pocket_type == "Fpocket" : 
            analysis.cart (dictionary_descriptor, list_without_ligand, pathDirectory.result (path_begin + "des_AutoSelect") +"des_select_desc") 
        else : 
            analysis.cart(dictionary_descriptor, list_with_ligand,  pathDirectory.result (path_begin + "des_AutoSelect") +"des_select_desc") 
            
    if radomForest : 
        
        path_begin = dataset_name + "/" + pocket_type + "/RANDOMFOREST/"
        analysis.radomForest(dictionary_descriptor, "global", pathDirectory.result (path_begin + "des_global") + "des_global")
        
    if regLog : 
        
        path_begin = dataset_name + "/" + pocket_type + "/RegLog/"
        analysis.glm(dictionary_descriptor, "global", pathDirectory.result (path_begin + "des_global") + "des_global")
        
    
    if correlation : 
        path_result = pathDirectory.result (dataset_name + "/" + pocket_type + "/CorDescriptor")
        analysis.correlationDescriptor (dictionary_descriptor, "global", path_result + "global_desc")
        
        # radi
        analysis.correlationDescriptor (dictionary_descriptor, "radi", path_result + "radi_desc")

        if "list_descriptor_ttest_signif" in locals () :
            analysis.correlationDescriptor (dictionary_descriptor, list_descriptor_ttest_signif,  path_result +"des_signiftest")
#      
        if pocket_type == "Fpocket" : 
            analysis.correlationDescriptor (dictionary_descriptor, list_without_ligand,  path_result +"des_select_desc") 
        else : 
            analysis.correlationDescriptor (dictionary_descriptor, list_with_ligand,  path_result +"des_select_desc") 
    
    if histogram : 
        path_begin = pathDirectory.result (dataset_name + "/" + pocket_type + "/Distribution")
        
        path_global_descriptor = writeFiles.globalDescriptors(dictionary_descriptor, path_begin + "global.data")
        runOtherProg.histogram (path_global_descriptor, path_begin, "drugg") 
        
     
    if SVM : 
        # path directory
        path_begin = pathDirectory.result (dataset_name + "/" + pocket_type + "/SVM")
        # global data for Rscript
        path_global_descriptor = writeFiles.globalDescriptors(dictionary_descriptor, path_begin + "global.data")
        
        
        if selection_variable == 1 :
            path_dir_parameter = pathDirectory.generatePath(path_begin + "selectedDesc/")
#             print path_dir_parameter
#            
            list_path_files = writeFiles.specificDescriptorbyData(dictionary_descriptor, "global", path_dir_parameter + "global_desc")
            list_descriptor_selected = parameters.retrieveDescriptorForSVMModel (dictionary_descriptor, path_dir_parameter, nb_descriptor = 3, debug = 0)   
            return
#        
        analysis.SVM (dictionary_descriptor, "global", pathDirectory.generatePath(path_begin + "global/"), "global", path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)

        
#       # descriptor selected
        if "list_descriptor_selected" in locals () :
            analysis.SVM (dictionary_descriptor, list_descriptor_selected, pathDirectory.generatePath(path_begin + "AutoSelected/"), "autoselected", path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
#        
#        # Descriptors imposed
#        if pocket_type == "Fpocket" : 
        analysis.SVM (dictionary_descriptor, list_without_ligand, pathDirectory.generatePath(path_begin + "WithoutLigand/"), "without_ligand", path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
#        else : 
        analysis.SVM (dictionary_descriptor, list_with_ligand, pathDirectory.generatePath(path_begin + "WithLigand/"), "with_ligand", path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
      
      
      
        if "list_descriptor_ttest_signif" in locals () :
            analysis.SVM (dictionary_descriptor, list_with_ligand, pathDirectory.generatePath(path_begin + "TtestSelected/"), "TtestSelected", path_file_global_descriptor = path_global_descriptor, name_dataset = dataset_name)
        
    #os.system ("rm " + path_file_color)
        

def analysisPocketEstimation (dictionnary_dataset, name_dataset, retrieve_type_pocket) :
    """
    Global analysis of pocket estimation (numbers pockets, ....)
    args: -> dictionary data set
          -> name data set to generate directory
    return: NONE, files results in directory results 
    
    """
    # number of pocket by PDB
    analysis.dataSetNumberPockets(dictionnary_dataset.keys (), name_dataset, retrieve_type_pocket)
    # distance between pockets
    checkPocket.distancePocket(dictionnary_dataset.keys(), pathDirectory.dataSet(name_dataset) ,pathDirectory.result(name_dataset + "/dataset/" + retrieve_type_pocket) + "distance_pocket")
    tool.generateFilePDBWithChain (dictionnary_dataset, pathDirectory.result(name_dataset + "/dataset/" + retrieve_type_pocket) + "list_PDB_chain", retrieve_type_pocket, name_dataset)
    tool.concatenePocket (dictionnary_dataset, pathDirectory.result(name_dataset + "/dataset/"+ retrieve_type_pocket) + "pocketsAll",name_dataset, retrieve_type_pocket, option_only_on_pocket = 1)
    


        
def  analysisDataset (dictionnary_dataset, name_dataset):
    """
    Analyse dataset DD between drug score and the confidence
    args: -> dictionnary dataset
          -> name dataset
    return: NULL draw histogram with drug score and confidence
    """
    analysis.resolution (dictionnary_dataset.keys(), name_dataset )
    if name_dataset == "DD" : 
        analysis.histCarac(dictionnary_dataset, name_dataset, "Drug Score")
        analysis.histCarac(dictionnary_dataset, name_dataset, "Confidence")
    
    
    
def calculDatasetDictionary (name_dataset, genrate_files = 0, debug = 0):
    
    path_file_krasowski = pathDirectory.dataSet() + "krasowski_2011_dataset.txt"
    path_file_Schmidtke = pathDirectory.dataSet() + "DD_Schmitke.csv"
    path_file_apo_Schmitke_clean = pathDirectory.dataSet() + "apoClean.txt"
    path_file_Perola = pathDirectory.dataSet() + "perola_phase3.txt"
    path_file_huang = pathDirectory.dataSet() + "apoCleanHuang.csv"
    dir_dataset = pathDirectory.dataSet ( name_dataset )
    
    if re.search ("krasowski", name_dataset) :
        dir_dataset = pathDirectory.dataSet ("krasowski")
        dico_dataset = dataSet.formatFileDataSetKrasowki2011(path_file_krasowski)
        if debug : print dico_dataset 
        if genrate_files : 
            for PDB_ID in dico_dataset.keys() :
                generateDataset (PDB_ID, dir_dataset)
            #identityCalcul.pasteFastaFileGlobal(dico_dataset.keys (), dir_dataset) # ancien protocole avec prog de leslie, maintenant generation direct -> water
            
        dataSet.compareDatasetKrasowskiSchmitke(dico_dataset, path_file_Schmidtke,dir_dataset, debug = 0)
        #if debug : print dico_dataset 
        dataSet.analysisLigandDataSet(dico_dataset, dir_dataset) # lipinski rule selected
        #if debug : print dico_dataset 
        dataSet.manualAnalysisForKrasowskiDataSet(dico_dataset)
        if debug : print len(dico_dataset.keys ()), ">TEST<"
        
        if name_dataset != "krasowski" :  # case retrieve only pocket train or test
            dataSet.diviseTrainTest (dico_dataset,name_dataset.split ("_")[-1] )
        
        
    elif name_dataset == "DD" : 
        dir_dataset = pathDirectory.dataSet ("DD")
        dico_dataset = dataSet.formatFileDataSchmitkeDD(path_file_Schmidtke)
        if genrate_files :
            for PDB_ID in dico_dataset.keys () : 
                generateDataset (PDB_ID, dir_dataset)
            identityCalcul.pasteFastaFileGlobal(dico_dataset.keys (), dir_dataset)
    
    elif name_dataset == "Perola" : 
        dir_dataset = pathDirectory.dataSet ("Perola")
        dico_dataset = dataSet.formatFileDataPerola(path_file_Perola)
        if genrate_files :
            for PDB_ID in dico_dataset.keys () : 
                generateDataset (PDB_ID, dir_dataset)
    
    
    elif name_dataset == "ApoForm" : 
        dir_dataset = pathDirectory.dataSet ("ApoForm")
        dico_dataset = dataSet.formatFileDataSchmitkeDD(path_file_Schmidtke)
        dico_dataset = dataSet.retrieveApoHoloFormDD (dico_dataset)
        if genrate_files :
            for PDB_ID in dico_dataset.keys () : 
                generateDataset (PDB_ID, dir_dataset)
            identityCalcul.pasteFastaFileGlobal(dico_dataset.keys (), dir_dataset)    
    
    elif name_dataset == "ApoForm128" : 
        dir_dataset = pathDirectory.dataSet ("ApoForm128")
        dico_dataset = dataSet.formatFileDataSchmitkeDD(path_file_Schmidtke)
        dico_dataset = dataSet.retrieveApoHoloFormDD (dico_dataset, l_PDB = dataSet.l_apo128)
        if genrate_files :
            for PDB_ID in dico_dataset.keys () : 
                generateDataset (PDB_ID, dir_dataset)
            identityCalcul.pasteFastaFileGlobal(dico_dataset.keys (), dir_dataset)  
    
    elif name_dataset == "ApoForm138" : 
        dir_dataset = pathDirectory.dataSet ("ApoForm138")
        dico_dataset = dataSet.formatFileDataSchmitkeDD(path_file_Schmidtke)
        dico_dataset = dataSet.retrieveApoHoloFormDD (dico_dataset, l_PDB = dataSet.l_apo138)
        if genrate_files :
            for PDB_ID in dico_dataset.keys () : 
                generateDataset (PDB_ID, dir_dataset)
            identityCalcul.pasteFastaFileGlobal(dico_dataset.keys (), dir_dataset)  
    
    
    elif name_dataset == "ApoFormClean" : 
        dir_dataset = pathDirectory.dataSet ("ApoFormClean")
        dico_dataset = dataSet.formatFileData(path_file_apo_Schmitke_clean)
        if genrate_files :
            for PDB_ID in dico_dataset.keys () : 
                generateDataset (PDB_ID, dir_dataset)
            identityCalcul.pasteFastaFileGlobal(dico_dataset.keys (), dir_dataset)       
    
    elif name_dataset == "ApoHuang" : 
        dir_dataset = pathDirectory.dataSet ("ApoHuang")
        dico_dataset = dataSet.formatFileDataSchmitkeDD(path_file_huang)
        dico_dataset = dataSet.retrieveApoHoloFormDD (dico_dataset)
        if genrate_files :
            for PDB_ID in dico_dataset.keys () : 
                generateDataset (PDB_ID, dir_dataset)
            identityCalcul.pasteFastaFileGlobal(dico_dataset.keys (), dir_dataset)     
          
    # Write dataset file (PDB with ligand pocket)       
    if debug : print dico_dataset 
#    writeFiles.datasetWithLigand(dico_dataset, dir_dataset + "dataSet.txt")
    return dico_dataset


def generateDataset (PDB_ID, dir_dataset):
    
    if downloadFile.importPDB( PDB_ID, dir_dataset) == 0 : 
        print "ERROR retrieve PDB files"
        return
    preparePDB.AllPDBSepareChain(PDB_ID, dir_dataset)
    #preparePDB.AllPDBProtonation(PDB_ID, dir_dataset) # -> protonation only protein.pdb not use now
    if downloadFile.importFasta(PDB_ID, dir_dataset) == 0 :
        #again
        downloadFile.importFasta(PDB_ID, dir_dataset)
        
    preparePDB.separeChainFasta (PDB_ID, dir_dataset)
    
    
def pocketEstimation (name_dataset, dictionary_dataset, pocket_type_retrieve, runFpocket = 1, file_dir_name = 1):
    """
    Estimation pocket
    args: -> name dataset
          -> dictionary dataset
          -> run Fpocket (run estimation)
          -> run Surflex (generate protomol)
          -> run NACCESS (generate protomol)
    return: NONE write files
    """
    
    dir_dataset = pathDirectory.dataSet(name_dataset)
    if pocket_type_retrieve == "proximity" :
        checkPocket.selectPocketProximity(dictionary_dataset, dir_dataset, name_dataset, pocket_type_retrieve, file_dir_name) 
    
    elif pocket_type_retrieve == "cavitator" : 
        runOtherProg.globalCavitator(dir_dataset, dictionary_dataset.keys(), debug = 1)
    
    elif pocket_type_retrieve == "DogSite" : 
        checkPocket.retrievePocketDogSite (dir_dataset, dictionary_dataset.keys(), debug = 1)
        checkPocket.selectPocketDogSite(dictionary_dataset, dir_dataset,  name_dataset, pocket_type_retrieve, file_dir_name = file_dir_name)
    
    else : 
        if runFpocket == 1 : 
            runOtherProg.globalFpocket(dir_dataset, dictionary_dataset.keys ())
        checkPocket.selectPocketFpocket(dictionary_dataset, dir_dataset,  name_dataset, pocket_type_retrieve, option_aggregation_pocket = 0, file_dir_name = file_dir_name)


def pocketEstimationApoForm (name_dataset, dictionary_dataset_apo, dictionary_dataset_holo, runFpocket = 1):
    
    # run Fpocket
    dir_dataset = pathDirectory.dataSet(name_dataset)
    if runFpocket == 1 : 
        runOtherProg.globalFpocket(dir_dataset, dictionary_dataset_apo.keys ())
    
    
    
    
    checkPocket.selectPocketFpocketApo(dictionary_dataset_apo, dictionary_dataset_holo, dir_dataset, name_dataset)


    
#     checkPocket.selectPocketFpocket(dictionary_dataset, dir_dataset,  name_dataset, pocket_type_retrieve, option_aggregation_pocket = 0, file_dir_name = file_dir_name)
#     
#     # superposed and run Fpocket one apo form
#     superposeStructure.superposeApoHolo(dictionary_dataset, name_dataset)
#     
#     # retrieve ligand in holo form
#     checkPocket.retrieveLigandFromHoloForm (dictionary_dataset, name_dataset)
#     
#     # select binding site
#     checkPocket.selectPocketFpocketApo(dictionary_dataset, name_dataset)
    
    

def generationProtomol (dictionary_dataset, pocket_type_retrieve, name_dataset):
    """
    Generate protomol with Surflex
    args: -> name data setname_dataset, pocket_type_retrieve
          -> dictionary data set
    return: -> generate Surflex protomol
    """
    dir_dataset = pathDirectory.dataSet(name_dataset)
    runOtherProg.protomolGenerationSurflexe (dictionary_dataset,  dir_dataset, pocket_type_retrieve, name_dataset )



def retrieveGlobalDescriptors (pocket_retrieve_type, protomol_type, dico_dataset, name_dataset, write_file = 1, calcul_descriptor = 0, option_separate = 1, compo_aa = 1):
    """
    Retrieve for every PDB in dico dataset values of descriptors
    args: -> type pocket retrieve
          -> protomol type
          -> dictionary dataset
          -> directory data set
    return: -> dictionary with descriptor
    """
    
    try : name_dataset = name_dataset.split ("_")[0]
    except  : pass
    
    
    
    dir_dataset = pathDirectory.dataSet(name_dataset)
    
    if calcul_descriptor == 1 : 
        getResidues.globalGetCompositionPocket (dico_dataset, dir_dataset, pocket_retrieve_type, name_dataset, debug = 1)
        descriptor.runDescriptorGlobal(dico_dataset, dir_dataset, protomol_type, pocket_retrieve_type, name_dataset, compo_aa = compo_aa)
    return loadDescriptors.generalLoadDescriptor(dico_dataset, protomol_type, pocket_retrieve_type, name_dataset, write_file=write_file, option_separate=option_separate)



def ACPDataset1BetwwenDataset2 (name_dataset1, name_dataset2, pocket_retrive_type1, protomol_type1, pocket_retrieve_type2, protomol_type2,  type_pocket = "all", list_descriptor = "global"):
    
    
    
    path_result = pathDirectory.result(name_dataset1 + "_" + name_dataset2)
    dataset1 = calculDatasetDictionary(name_dataset1, 0)
    dataset2 = calculDatasetDictionary(name_dataset2, 0)
    
    dico_data1= retrieveGlobalDescriptors (pocket_retrive_type1, protomol_type1, dataset1, name_dataset1, write_file = 0, calcul_descriptor = 0, option_separate = 1)
    dico_data2= retrieveGlobalDescriptors (pocket_retrieve_type2, protomol_type2, dataset2, name_dataset2, write_file = 0, calcul_descriptor = 0, option_separate = 1)
    
#    print dico_data1["No-Druggable"]["PDB"]
#    print dico_data2["No-Druggable"]["PDB"]
#    list_union = list(set(dico_data2["No-Druggable"]["PDB"])& set( dico_data1["No-Druggable"]["PDB"]))
#    print list_union
    
    path_file_data1 = path_result + name_dataset1 + "1"
    path_file_data2 = path_result + name_dataset2 + "2"
    path_file_result = path_result + name_dataset1 + "_" + name_dataset2
    path_file_color = writeFiles.colorACPFile (dico_data1)
    
    if list_descriptor == "global" : 
        writeFiles.globalDescriptors(dico_data1, path_file_data1, type_pocket = type_pocket)
        writeFiles.globalDescriptors(dico_data2, path_file_data2, type_pocket = type_pocket)
    else : 
        writeFiles.specificDescriptor(dico_data1, list_descriptor, path_file_data1, debug = 1 )
        writeFiles.specificDescriptor(dico_data2, list_descriptor, path_file_data2, debug = 1 )
    
    runOtherProg.ACPDataset (path_file_data1, path_file_data2, path_file_result)


def MDSDataset1BetwwenDataset2 (name_dataset1, name_dataset2, pocket_retrive_type1, protomol_type1, pocket_retrieve_type2, protomol_type2,  type_pocket = "all", list_descriptor = "global"):
    
    
    
    path_result = pathDirectory.result(name_dataset1 + "_" + name_dataset2)
    dataset1 = calculDatasetDictionary(name_dataset1, 0)
    dataset2 = calculDatasetDictionary(name_dataset2, 0)
    
    dico_data1= retrieveGlobalDescriptors (pocket_retrive_type1, protomol_type1, dataset1, name_dataset1, write_file = 0, calcul_descriptor = 0, option_separate = 1)
    dico_data2= retrieveGlobalDescriptors (pocket_retrieve_type2, protomol_type2, dataset2, name_dataset2, write_file = 0, calcul_descriptor = 0, option_separate = 1)
    
#    print dico_data1["No-Druggable"]["PDB"]
#    print dico_data2["No-Druggable"]["PDB"]
#    list_union = list(set(dico_data2["No-Druggable"]["PDB"])& set( dico_data1["No-Druggable"]["PDB"]))
#    print list_union
    
    path_file_data1 = path_result + name_dataset1 + "1"
    path_file_data2 = path_result + name_dataset2 + "2"
    path_file_result = path_result + name_dataset1 + "_" + name_dataset2
    path_file_color = writeFiles.colorACPFile (dico_data1)
    
    if list_descriptor == "global" : 
        writeFiles.globalDescriptors(dico_data1, path_file_data1, type_pocket = type_pocket)
        writeFiles.globalDescriptors(dico_data2, path_file_data2, type_pocket = type_pocket)
    else : 
        writeFiles.specificDescriptor(dico_data1, list_descriptor, path_file_data1, debug = 1 )
        writeFiles.specificDescriptor(dico_data2, list_descriptor, path_file_data2, debug = 1 )
    
    runOtherProg.MDSDataset (path_file_data1, path_file_data2, path_file_result)









def correlationDescriptorsByTypePocket (name_dataset, pocket_type1, protomol_type1, pocket_type2, protomol_type2) : 
    """Correlation between two data of descriptor"""
    
    dico_dataset = calculDatasetDictionary(name_dataset, 0)
    path_dir_result = pathDirectory.result(name_dataset + "/" + pocket_type1 + "VS" + pocket_type2)
    
    dico_descriptor1 = retrieveGlobalDescriptors (pocket_type1, protomol_type1, dico_dataset, name_dataset)
    dico_descriptor2 = retrieveGlobalDescriptors (pocket_type2, protomol_type2, dico_dataset, name_dataset)
    
    writeFiles.globalDescriptors(dico_descriptor1, path_dir_result + "desc1")
    writeFiles.globalDescriptors(dico_descriptor2, path_dir_result + "desc2")
    
    path_correlation = analysis.correlationDotChart (dico_descriptor1, dico_descriptor2, path_dir_result)
    runOtherProg.dotChart(path_correlation)
    runOtherProg.plotCorrelation (path_dir_result + "desc1", path_dir_result + "desc2", path_dir_result + "correlationbydescriptor")


def validation (name_dataset_train, name_dataset_test, pocket_retrieve_type, protomol_type, pvalue = 0.1, cor_value = 1):
    """
    A revoir car plus utilise, moyen de passer avec un save.image -----> a voir
    
    
    """
    # file color for ACP and path
    path_result = pathDirectory.result ("modelDruggability" + "/" + name_dataset_train + "_" + name_dataset_test)

    # load dataset
    dico_dataset_train = calculDatasetDictionary(name_dataset_train, 0)
    dico_dataset_test = calculDatasetDictionary(name_dataset_test, 0)
    
    # load descriptor
    dico_descriptors_train = retrieveGlobalDescriptors (pocket_retrieve_type, protomol_type, dico_dataset_train, name_dataset_train, write_file = 0, calcul_descriptor = 0, option_separate = 1)
    dico_descriptors_test = retrieveGlobalDescriptors (pocket_retrieve_type, protomol_type, dico_dataset_test, name_dataset_test, write_file = 0, calcul_descriptor = 0, option_separate = 1)
    
    # color ACP
    writeFiles.colorACPFile (dico_descriptors_train) # color for ACP
    
    # global file for R script
    path_file_global_train = writeFiles.globalDescriptors(dico_descriptors_train, "global.data") # for Ttest and histogram, open in R script
    path_file_global_test = writeFiles.globalDescriptors(dico_descriptors_test, "global_test.data")
    
    # retrieve descriptor significant mean score
    if len (dico_descriptors_train["Druggable"].keys()) > 1 and len (dico_descriptors_train["No-Druggable"].keys()) > 1 :  
        dico_ttest = analysis.ttest(dico_descriptors_train["Druggable"], dico_descriptors_train["No-Druggable"], path_result + "ttest_" + name_dataset_train )
        # descriptor positive ttest between druggable and non-druggable
        dico_IC_signif = analysis.retrieveDescriptorSignifTtest (dico_ttest, pvalue = pvalue)

    list_descriptor_ttest_signif = dico_IC_signif.keys ()
       
    # write data with significatives descriptor
#    path_signif_descriptor_train = writeFiles.specificDescriptor(dico_descriptors_train, list_descriptor_ttest_signif, path_result + name_dataset_train+ "_train" )
#    path_signif_descriptor_test = writeFiles.specificDescriptor(dico_descriptors_test, list_descriptor_ttest_signif, path_result + name_dataset_test+ "_test" )
    
    # a voir 
    
    
    # parameter in LDA value correlation -> elimcor
    #runOtherProg.accTrainTest(path_signif_descriptor_train, path_signif_descriptor_test, 0.4,  0.02)
    #return
    list_descriptor_selected, IC_descriptor = analysis.LDAGlobal (dico_descriptors_train, list_descriptor_ttest_signif, path_result, "model", name_dataset = name_dataset_train)
    
    # descriptor signif betwenn good predicted and bad predicted
    path_descriptor = path_result + "descriptor_signif_good_bad_predicted"
    filout_descriptor =open (path_descriptor, "w")
    for descriptor in IC_descriptor.keys () : 
        if descriptor == "%_ATOM_CONVEXE" : 
            descriptor = "X._ATOM_CONVEXE"
        elif descriptor == "Mean_alpha-sphere_SA" : 
            descriptor = "Mean_alpha.sphere_SA"
        filout_descriptor.write (descriptor + "\n")
    filout_descriptor.close ()
    
    dico_test_corecting = tool.selectData (dico_descriptors_test, IC_descriptor)
    
    
    list_predicted = analysis.LDAmodel(dico_descriptors_train, dico_descriptors_test, dico_test_corecting, list_descriptor_selected, path_file_global_train, path_file_global_test, 1, path_result, path_descriptor)
    
#    os.system ("rm " + path_file_global_test)
#    os.system ("rm " + path_file_global_train)
    
#    analysis.scoreBadPredict (list_predicted, dico_dataset_train, dico_dataset_test)


def applyModel(path_file_model, path_file_descriptor,  path_file_result, family_barplot = 0, name_dataset = "") :
    """
    Apply model LDA one new data
    args: -> path file model
          -> path file descriptors (new data)
          -> path file result
    return: path file result
    """
    if os.path.exists(path_file_result) : 
        os.remove(path_file_result)
    
    path_proba, path_result = runOtherProg.predictLDA (path_file_model, path_file_descriptor, path_file_result)
    
    if family_barplot : 
        dico_data = calculDatasetDictionary(name_dataset,0)
        filin = open (path_proba,"r")
        l_PDBscore = filin.readlines()[1:]
        filin.close ()
        filout = open (path_proba, "w")
        for PDBscore in l_PDBscore : 
            PDB = PDBscore.split (" ")[0].replace("\"", "")
            filout.write(str(PDB) + "\t" + str(PDBscore.split (" ")[1]) + "\t" + str(dico_data[PDB]["Protein name"]) + "\t" + str(dico_data[PDB]["druggability"]) + "\n")
        filout.close()
        
        runOtherProg.boxplotFamily (path_proba)
    
    
    return path_file_result
    
    
    



def correlationDescriptor (name_dataset, list_desc, pocket_retrive_type, protomol_type):
    
    dico_dataset = calculDatasetDictionary(name_dataset, 0)
    dico_descriptors = retrieveGlobalDescriptors (pocket_retrive_type, protomol_type, dico_dataset, name_dataset, write_file = 0, calcul_descriptor = 0, option_separate = 1)
    
    
    analysis.correlationBetween2Descriptor(dico_descriptors, list_desc)
    
    

def overlapTwoPockets (name_dataset, pocket_retrieve_type_1, pocket_retrieve_type_2): 
    
    
    path_dir_result = pathDirectory.result("overlapPocket" + str (pocket_retrieve_type_1) + "_" + str (pocket_retrieve_type_2))
    path_filout = path_dir_result + "histogramScore"
    filout = open(path_filout, "w")
    filout.write ("PDB\tPO\tMO\tRO\tcomon\tnb1\tnb2\n")
    dico_dataset = calculDatasetDictionary(name_dataset, 0)
    dico_value = {}
    
    
    for PDB_ID in dico_dataset.keys () : 
        dico_value[PDB_ID] = []
        list_dir_pocket1 = pathDirectory.generateListDirPocket (PDB_ID, pocket_retrieve_type_1, name_dataset)
        list_dir_pocket2 = pathDirectory.generateListDirPocket (PDB_ID, pocket_retrieve_type_2, name_dataset)
        if not list_dir_pocket1  or not  list_dir_pocket2  : 
            continue
        
        
        for dir_pocket1 in list_dir_pocket1 : 
            for dir_pocket2 in list_dir_pocket2 : 
                path_file_pocket1 = pathDirectory.searchPocketAtomFpocket(dir_pocket1)
                path_file_pocket2 = pathDirectory.searchPocketAtomFpocket(dir_pocket2)
                
                p_file_pocket1_asa = pathDirectory.searchPocketACC (dir_pocket1)
                p_file_pocket2_asa = pathDirectory.searchPocketACC (dir_pocket2)
                
                score_overlap = overlapPocket.scoreOverlap(path_file_pocket1, path_file_pocket2)
                nb_commun_atom, nb_atom1, nb_atom2 = overlapPocket.communAtom(path_file_pocket1, path_file_pocket2)
                score_MO = overlapPocket.MO(path_file_pocket1, path_file_pocket2, p_file_pocket1_asa, p_file_pocket2_asa)
                score_RO = overlapPocket.MO(path_file_pocket2, path_file_pocket1, p_file_pocket2_asa, p_file_pocket1_asa)
                dico_value[PDB_ID].append(score_overlap)
                dico_value[PDB_ID].append(score_MO)
                dico_value[PDB_ID].append(score_RO)
                dico_value[PDB_ID].append(nb_commun_atom)
                dico_value[PDB_ID].append(nb_atom1)
                dico_value[PDB_ID].append(nb_atom2)
        
                   
        filout.write (str(PDB_ID) + "\t" + str (dico_value[PDB_ID][0]) + "\t" + str (dico_value[PDB_ID][1])+ "\t" + str (dico_value[PDB_ID][2]) + "\t" + str (dico_value[PDB_ID][3]) + "\t" + str (dico_value[PDB_ID][4]) + "\t" + str (dico_value[PDB_ID][5]) + "\n")
    
   
    filout.close () 
    runOtherProg.histScoreByPDB (path_filout)
        
    

def retrieveAccessibilityPocket(name_dataset, pocket_type) :
    
    path_dir_result = pathDirectory.result("accesibility" + str (pocket_type))
    path_filout = path_dir_result + "accessibility_pocket.out"
    filout = open(path_filout, "w")
    dico_dataset = calculDatasetDictionary(name_dataset, 0)
    
    for PDB_ID in dico_dataset.keys () : 
        path_pocket_asa = pathDirectory.searchPocketAtomASA(PDB_ID, name_dataset, pocket_type)
        if path_pocket_asa == None : continue
        acc = descriptorEnergy.accessibilityPocketbyAtom(path_pocket_asa)
        filout.write (PDB_ID + "\t" + str (acc) + "\n")
    
    filout.close ()
        
        
        
        
        

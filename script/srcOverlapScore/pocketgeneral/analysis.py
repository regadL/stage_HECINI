"""
BORREL Alexandre
04-2012
"""
import pathDirectory
import runOtherProg
import tool
import writeFiles
import downloadFile
import parseFpocket
import parameters
from PDB import *
import globalFonction
import parsePDB


import os
from re import search
from numpy import arange, mean
from copy import deepcopy
from math import sqrt
import loadDescriptors
import dataSet
import superposeStructure


def ttest (dico_descriptorIn1, dico_descriptorIn2, path_filout,  type_ttest = "drug" ,sorting = 0, option_cor = 0, debug = 0):
    """
    Calcul Ttest with different descriptors between 2 calcul type
    args:-> dictionary 1 with descriptor
         -> dictionary 2 with descriptor
         -> path file out
         -> option correlation
    return: -> dictionary with major key -> descriptor andvalue test
    """
    dico_out = {}
    
    # copy dictionary because modify dictionary
    dico_descriptor1 = deepcopy(dico_descriptorIn1)
    dico_descriptor2 = deepcopy(dico_descriptorIn2)
    
    try : 
        del dico_descriptor1["PDB"]
        del dico_descriptor1["data"]
    except :
        pass
    
    try :
        del dico_descriptor2["PDB"]
        del dico_descriptor2["data"]
    except :
        pass
    
    # check if last dictionary keys
    if "Druggable" in dico_descriptor1.keys () and "Druggable" in dico_descriptor2.keys () :
        dico_descriptor1 = tool.fusionDictionaryPocketType (dico_descriptor1, data_separated = 0)
        del dico_descriptor1["PDB"]
        dico_descriptor2 = tool.fusionDictionaryPocketType (dico_descriptor2, data_separated = 0)
        del dico_descriptor2["PDB"]
    
    for type_descriptor in dico_descriptor1.keys () : 
        for descriptor in dico_descriptor1[type_descriptor].keys () :
            dico_out[descriptor] = {} 
            path_tem1 = path_filout + "_temp1"
            path_tem2 = path_filout + "_temp2"
            writeFiles.writeListValues(dico_descriptor1[type_descriptor][descriptor], path_tem1)
            writeFiles.writeListValues(dico_descriptor2[type_descriptor][descriptor], path_tem2)
            
            if "NA" in dico_descriptor1[type_descriptor][descriptor] or "-nan" in dico_descriptor1[type_descriptor][descriptor]: 
                tool.suppNA(dico_descriptor1[type_descriptor][descriptor])
            if "NA" in dico_descriptor2[type_descriptor][descriptor] or "-nan" in dico_descriptor2[type_descriptor][descriptor]: 
                tool.suppNA(dico_descriptor2[type_descriptor][descriptor])
            
            # number value   
            nb_value = len (dico_descriptor1[type_descriptor][descriptor]) + len(dico_descriptor2[type_descriptor][descriptor])
            if debug :
                print type_descriptor, descriptor
                print dico_descriptor1[type_descriptor][descriptor]
                print dico_descriptor2[type_descriptor][descriptor]
                
            try : descriptor = descriptor.replace( " ", "_" )
            except : pass
            dico_out[descriptor] = runOtherProg.runRscriptTtest( path_tem1, path_tem2 )
            dico_out[descriptor]["nb values"] = int(nb_value)
            dico_out[descriptor]["nb drugg"] = int(len (dico_descriptor1[type_descriptor][descriptor]))
            dico_out[descriptor]["nb no drugg"] = int(len(dico_descriptor2[type_descriptor][descriptor]))
    
    # rm path tem
    os.system ("rm " + path_tem1)
    os.system ("rm " + path_tem2)
    if debug : print dico_out
    writeFiles.resultTtest(dico_out, path_filout, type_ttest, sorting, option_cor) 
    return dico_out

            
def correlationVolumeVSVolumeFpocket (dico_descriptors, list_PDB, name_dataSet, pocket_retrieve_type) : 
    """
    Correlation volume
    args: -> dictionary descriptors
          -> list PDB
          -> path filout
    return: NONE
    Draw plot with volume Fpocket vs Surflexe
    """
    path_filout = pathDirectory.result(name_dataSet + "/" + pocket_retrieve_type) + "cor_volume"
    filout = open (path_filout, "w")
    filout.write ("PDB\tVolume Surflexe\tVolume Fpocket\n") # header for R script
    list_type_pocket = dico_descriptors.keys ()
    list_type_pocket.remove ("data")
    for type_pocket in list_type_pocket:
        for i in xrange (0,len(dico_descriptors[type_pocket]["volume"]["volume"])) : 
            try : filout.write (dico_descriptors[type_pocket]["PDB"][i] + "\t" + dico_descriptors[type_pocket]["volume"]["volume"][i] + "\t" + dico_descriptors[type_pocket]["fpocket"]["Real volume"][i] + "\n")
            except : filout.write (dico_descriptors[type_pocket]["PDB"][i] + "\t" + dico_descriptors[type_pocket]["volume"]["volume"][i] + "\t" + dico_descriptors[type_pocket]["fpocket"]["pock_vol"][i] + "\n")
    filout.close () 
    
    runOtherProg.runRcorrelation(path_filout)
    os.system ("rm " + path_filout)
    

def dataSetNumberPockets (list_PDB, name_dataset, retrieve_type_pocket):
    """
    Retrieve by protein the number of pocket
    args: -> list PDB
          -> name data set
    return: NONE write file pocket analysis
    """
    
    path_dir_descriptor = pathDirectory.descriptor(name_dataset + "/" + retrieve_type_pocket)
    
    filout = open (pathDirectory.result(name_dataset + "/dataset/" + retrieve_type_pocket) + "pocket_analysis.txt", "w")
    for PDB_ID in list_PDB : 
        nb_pocket = 0
        path_result_PDB = path_dir_descriptor + PDB_ID
        list_files = os.listdir(path_result_PDB)
        for file_PDB in list_files : 
            if search ("^pocket",file_PDB) : 
                nb_pocket = nb_pocket + 1
        filout.write (str (PDB_ID) + "\t" + str (nb_pocket) + "\n")
    filout.close ()
        

def specificACP (descriptor_in , dico_descriptor, path_filout, mainACP = "Title"):
    
    if descriptor_in == "global" : 
        writeFiles.globalDescriptors(dico_descriptor,path_filout )
        runOtherProg.runRACP(path_filout, mainACP)
        
    else : 
        writeFiles.specificDescriptor(dico_descriptor, descriptor_in, path_filout )
        runOtherProg.runRACP(path_filout, mainACP)
    
    #os.system ("rm " + path_filout)
       



def correlationDotChart (dico_descriptor1, dico_descriptor2, path_directory, debug = 0):
    
    path_filout = path_directory + "correlation.data"
    filout = open (path_filout, "w")
    
    list_cor_write = ""
    l_pvalue = ""
    list_descriptor_write = ""
    
    # check if last dictionary keys
    if "Druggable" in dico_descriptor1.keys () and "No-Druggable" in dico_descriptor1.keys () :
        dico_descriptor1 = tool.fusionDictionaryPocketType (dico_descriptor1)
    elif "data" in dico_descriptor1.keys () : 
        del dico_descriptor1["data"]
        
    if "Druggable" in dico_descriptor2.keys () and "No-Druggable" in dico_descriptor2.keys () : 
        dico_descriptor2 = tool.fusionDictionaryPocketType (dico_descriptor2)
    elif "data" in dico_descriptor2.keys () :
        del dico_descriptor1["data"]
    
    
    tool.balanceSamePDB (dico_descriptor1,dico_descriptor2)
    
    for type_descriptor in dico_descriptor1.keys () : 
        if type_descriptor == "data" : continue
        for descriptor in dico_descriptor1[type_descriptor].keys () :
            path_tem1 = path_directory + "temp1"
            path_tem2 = path_directory + "temp2"
            try : 
                writeFiles.writeListValues(dico_descriptor1[type_descriptor][descriptor], path_tem1)
                writeFiles.writeListValues(dico_descriptor2[type_descriptor][descriptor], path_tem2)
                list_descriptor_write = list_descriptor_write + str(descriptor) + "\t"
                
                #print len (dico_descriptor1[type_descriptor][descriptor]), len (dico_descriptor2[type_descriptor][descriptor]),
            except : 
                continue
            
            if "NA" in dico_descriptor1[type_descriptor][descriptor] : 
                tool.suppNA(dico_descriptor1[type_descriptor][descriptor])
            if "NA" in dico_descriptor2[type_descriptor][descriptor] : 
                tool.suppNA(dico_descriptor2[type_descriptor][descriptor])
                
            if debug :
                print type_descriptor, descriptor
                print dico_descriptor1[type_descriptor][descriptor]
                print dico_descriptor2[type_descriptor][descriptor]
                
            dico_ttest = runOtherProg.runRscriptTtest( path_tem1, path_tem2 )
            list_cor_write = list_cor_write + str (dico_ttest["corr"]) + "\t"
            l_pvalue = l_pvalue + str (dico_ttest["p-value_cor"]) + "\t"
            
            
    filout.write (list_descriptor_write[:-1] + "\n")
    filout.write (list_cor_write [:-1]+ "\n")
    filout.write (l_pvalue [:-1]+ "\n")
    filout.close ()
    return path_filout
    
    


def fpocketParameters (list_PDB, step = 0.5, run = 1, debug = 1):
    
    dir_file_PDB = pathDirectory.FpocketTest()
    if run : 
        downloadFile.importPDB(list_PDB, dir_file_PDB, dir_by_PDB=0, debug = 1)
    
        for m in arange(2.5,3.6, step) :
            m = str (m)
            if m[-2:] == ".0" :
                m = m.split (".")[0]
            filout_nb_pocket = open (dir_file_PDB + "analysisNbPocket_" + str (m), "w")
            filout_vol = open (dir_file_PDB + "analysisVol_" + str (m), "w")
            for M in arange(5.5,6.6, step) :
                print M
                filout_nb_pocket.write (str (M) + "\t")
                filout_vol.write (str (M) + "\t")
                list_vol = []
                list_nb_pocket = []
                for PDB_ID in list_PDB : 
                    runOtherProg.runFpocket(dir_file_PDB + PDB_ID.upper () + ".pdb", m = m, M = M)
                    dir_out_fpocket = dir_file_PDB + PDB_ID.upper () + "_out/"
                    parseFpocket.volumeNbPocket (list_vol, list_nb_pocket, dir_out_fpocket)
                    os.system ("rm -r " + dir_file_PDB + "*out")
                # mean on list PDB
                filout_nb_pocket.write (str (mean(list_nb_pocket)) + "\n")
                filout_vol.write (str (mean(list_vol)) + "\n")
            filout_nb_pocket.close ()
            filout_vol.close ()
    
    runOtherProg.runRhistoFpocket (dir_file_PDB + "analysisNbPocket", "pocket", begin = 2.5, end = 3.5, step = step)
    runOtherProg.runRhistoFpocket (dir_file_PDB + "analysisVol", "Volume", begin = 2.5, end = 3.5, step = step)
            
            
            

def LDAGlobal (dico_descriptors, list_descriptors, path_directory, file_out,  name_barplot = "0", test_train_data = 1, graph = 1, bad_analysis = 1, name_dataset = "", impose_nb_descriptor = 0, path_file_global_descriptor = ""):
    
    path_filout = path_directory + file_out

    ###############
    # write files #
    ###############
    # LOO
    
    if list_descriptors == "global" :
        path_file_every_pocket =writeFiles.globalDescriptors(dico_descriptors, path_filout + "_loo")
    else : 
        path_file_every_pocket = writeFiles.specificDescriptor(dico_descriptors, list_descriptors, path_filout + "_loo" )
    
    # train and test
    if test_train_data  :
        path_files_train_test = writeFiles.specificDescriptorbyData (dico_descriptors, list_descriptors, path_filout)
    
    ###############
    #   RUN LDA   #
    ###############
    # run LDA loo  -> elimcor -> 0 because change selected variables      
    
    runOtherProg.bartlett (path_filin = path_file_every_pocket, path_filout = path_filout + "_bartlett_validity") 
    runOtherProg.lda (path_filin1 = path_file_every_pocket, path_filout = path_filout, path_file_global = path_file_global_descriptor, leave_one_out = 1, name_barplot = name_barplot, elimcor = 0, graph=graph)
    
    # run LDA train and test
    if test_train_data == 1 : 
        if not os.path.exists(path_filout + "_train") and not os.path.exists(path_filout + "_test"): # case elimcor before
            writeFiles.specificDescriptorbyData (dico_descriptors, list_descriptors, path_filout)
        runOtherProg.lda (path_filin1 =  path_files_train_test [0], path_filin2 =  path_files_train_test [1], path_filout = path_filout, path_file_global = path_file_global_descriptor, leave_one_out = 0, name_barplot = name_barplot, elimcor = 0, graph=graph)
   
    # curve ROC
    runOtherProg.ldaROCCurve (path_filin1 =  path_files_train_test [0], path_filin2 =  path_files_train_test [1], path_filout = path_filout)
   
    if bad_analysis : 
        dir_bad_predict = pathDirectory.generatePath(path_directory + "BadPredict/")
        dico_ttest = ttestBetwwenBadGoodPredict ("good_predict", "bad_predict", dico_descriptors, dir_bad_predict, "global_good_bad", rm = 0)
        ttestBetwwenBadGoodPredict ("good_predict_no_drug", "bad_predict_no_drug", dico_descriptors, dir_bad_predict, "no_drug")
        ttestBetwwenBadGoodPredict ("good_predict_drug", "bad_predict_drug", dico_descriptors, dir_bad_predict, "drug")
        badPrediction ("good_predict", "bad_predict", dico_descriptors, dir_bad_predict)
        dictionary_IC = retrieveDescriptorSignifTtest (dico_ttest, pvalue = 0.005)

    
    list_file_barplot = pathDirectory.searchModel (path_directory)
    
    # a modifier construire une seule structure avec direct dataset et descripteur
#    scoreBadPredict (list_file_barplot,  globalFonction.calculDatasetDictionary(name_dataset, 0), globalFonction.calculDatasetDictionary(name_dataset, 0))
    return list_descriptors, dictionary_IC


def badPrediction (path_file_good, path_file_bad, dico_descriptor, path_dir_bad_predict) : 

    list_PDB_good = tool.loadFileWithList(path_file_good)
    list_PDB_bad = tool.loadFileWithList(path_file_bad)
    
    dico_bad = loadDescriptors.loadDescriptorSpecificPDBID(dico_descriptor, list_PDB_bad)
    path_file_bad_global = writeFiles.globalDescriptors(dico_bad, path_dir_bad_predict + "bad_predicted" )

    # run to draw tree plot by cart    
    runOtherProg.cartPlot (path_file_bad_global, path_dir_bad_predict + "cart_bad")
    
    os.system ("rm " + path_file_good)
    os.system ("rm " + path_file_bad)
    os.system ("rm " + path_file_bad_global)



def ldaLeaveOneOut (dico_descriptors, list_descriptor, path_filin, elimcor = 0) : 

    if list_descriptor == "global" :
        writeFiles.globalDescriptors(dico_descriptors,path_filin )
    else :
        writeFiles.specificDescriptor(dico_descriptors, list_descriptor, path_filin )
    
    if elimcor : 
        elimcor = parameters.onlyACCLoo(path_filin)
        os.system ("rm " + path_filin + "*.png")
        
    runOtherProg.lda (path_filin1 = path_filin, path_filout = path_filin + "_loo" , leave_one_out = 1, elimcor = elimcor)

def ldaData (dico_descriptors, list_descriptor, path_filout) : 
    
    writeFiles.specificDescriptorbyData (dico_descriptors, list_descriptor, path_filout )
    runOtherProg.lda (path_filin1 = path_filout + "_train", path_filin2 = path_filout + "_test", leave_one_out = 0)
    

def cart(dico_descriptors, list_descriptor, path_filout, elimcor=0, train_test_data = 1) : 
    """
    Run cart
    args: -> dictionary descriptors
          -> list descriptors
          -> path file out
    Return: NONE (run R cart)
    """
    
    #print path_filout
    if list_descriptor == "global" : 
        path_file_global = writeFiles.globalDescriptors(dico_descriptors,path_filout )
    else : 
        path_file_global = writeFiles.specificDescriptor(dico_descriptors, list_descriptor, path_filout )
    
    if train_test_data == 1 : 
        path_file_by_data = writeFiles.specificDescriptorbyData (dico_descriptors, list_descriptor, path_filout)
        runOtherProg.Rcart (path_file_global, 0, 0, elimcor) # path create in R script !!!! a refaire pour homogeniser
        runOtherProg.Rcart (path_file_global, path_file_by_data[0],  path_file_by_data[1], elimcor)
        
    else : 
        runOtherProg.Rcart (path_file_global, 0, 0, elimcor)
    #os.system ("rm " + path_filout + "_train " + path_filout + " " + path_filout + "_test")


def radomForest (dico_descriptors, list_descriptor, path_filout, elimcor=0, train_test_data = 1) : 
    """
    Run RF
    args: -> dictionary descriptors
          -> list descriptors
          -> path file out
    Return: NONE (run R RF)
    """
    
    #print path_filout
    if list_descriptor == "global" : 
        path_file_global = writeFiles.globalDescriptors(dico_descriptors,path_filout )
    else : 
        path_file_global = writeFiles.specificDescriptor(dico_descriptors, list_descriptor, path_filout )
    
    if train_test_data == 1 : 
        path_file_by_data = writeFiles.specificDescriptorbyData (dico_descriptors, list_descriptor, path_filout)
        runOtherProg.RRandomForest (path_file_global, 0, 0, elimcor) # path create in R script !!!! a refaire pour homogeniser
        runOtherProg.RRandomForest (path_file_global, path_file_by_data[0],  path_file_by_data[1], elimcor)
        
    else : 
        runOtherProg.RRandomForest (path_file_global, 0, 0, elimcor)
    #os.system ("rm " + path_filout + "_train " + path_filout + " " + path_filout + "_test")



def glm(dico_descriptors, list_descriptor, path_filout, elimcor=0, train_test_data = 1) : 
    
    #print path_filout
    if list_descriptor == "global" : 
        path_file_global = writeFiles.globalDescriptors(dico_descriptors,path_filout )
    else : 
        path_file_global = writeFiles.specificDescriptor(dico_descriptors, list_descriptor, path_filout )
    
    if train_test_data == 1 : 
        path_file_by_data = writeFiles.specificDescriptorbyData (dico_descriptors, list_descriptor, path_filout)
        runOtherProg.Rglm (path_file_global, 0, 0, elimcor) # path create in R script !!!! a refaire pour homogeniser
        runOtherProg.Rglm (path_file_global, path_file_by_data[0],  path_file_by_data[1], elimcor)
        
    else : 
        runOtherProg.Rglm (path_file_global, 0, 0, elimcor)

    

def correlationDescriptor (dico_descriptors, list_descriptor, path_filout, nb_color = 6) : 
    """
    correlation betwen descriptors
    args: -> dictionary with descriptors
          -> list descriptor
          -> path file out
          -> nb color in card
    return: NONE (run R card script)
    """
    
    if list_descriptor == "global" : 
        writeFiles.globalDescriptors(dico_descriptors,path_filout )
    else : 
        writeFiles.specificDescriptor(dico_descriptors, list_descriptor, path_filout )
        
    runOtherProg.MDS (path_filout)  
    runOtherProg.matrixCorr (path_filout, nb_color)    
    
#     os.system ("rm " + path_filout)
    
    
def retrieveDescriptorSignifTtest (dico_ttest, pvalue= 0.05, debug = 0) : 
    """
    Retrieve with ttest dictionary the significative descriptor
    args: -> dictionary ttest
          -> option selected value
          -> pvalues
    return: IC for interest descriptor and write list descriptors
    """
    
    dico_IC = {}
    for descriptor in dico_ttest.keys () :
        if dico_ttest[descriptor]["p-value"] < pvalue : 
            dico_IC[descriptor] = {}
            dico_IC[descriptor]["Borne inf"] = dico_ttest[descriptor]["no_drugg"] - 1.96 * (dico_ttest[descriptor]["sd_no_drugg"] / sqrt(dico_ttest[descriptor]["nb no drugg"]))
            dico_IC[descriptor]["Borne sup"] = dico_ttest[descriptor]["no_drugg"] + 1.96 * (dico_ttest[descriptor]["sd_no_drugg"] / sqrt(dico_ttest[descriptor]["nb no drugg"]))
    
    if debug : print dico_IC
    return dico_IC
        


def ttestBetwwenBadGoodPredict (path_file_good, path_file_bad, dico_descriptor, path_dirout, type_data, rm = 1) :
    
    list_PDB_good = tool.loadFileWithList(path_file_good)
    list_PDB_bad = tool.loadFileWithList(path_file_bad)
    
    dico_good = loadDescriptors.loadDescriptorSpecificPDBID(dico_descriptor, list_PDB_good)
    dico_bad = loadDescriptors.loadDescriptorSpecificPDBID(dico_descriptor, list_PDB_bad)
    
    out_ttest = ttest (dico_good, dico_bad, path_dirout + type_data + ".result", type_ttest = "predicted", sorting = 1, option_cor = 0, debug = 0)

    if rm : 
        os.system ("rm " + path_file_good)
        os.system ("rm " + path_file_bad)
        
    return out_ttest

def resolution (list_PDB, name_dataset) : 
    """
    Analysis resolution for every PDB
    arg: -> list of PDB ID
         -> name dataset for generate path dataset
    return: NONE (draw histogram)
    """
    path_filout = pathDirectory.result(name_dataset + "/dataset" ) + "resolution" 
    filout = open (path_filout, "w")
    for PDB_ID in list_PDB :
        path_PDB = pathDirectory.dataSet(name_dataset) + PDB_ID + "/" + PDB_ID + ".pdb"
        try : 
            PDB_parsed = PDB (path_PDB)
            resol = PDB_parsed.resolution(verbose = 0)
            filout.write (str (resol) + "\n")
        except:
            pass
    filout.close ()
    
    runOtherProg.runRscriptHisto(path_filout, "Resolution")
    
def histCarac (dico_dataset, name_dataset, in_keys) : 
    """
    Analysis drug score for every PDB in DD
    arg: -> dictionary dataset
         -> name dataset for generate path dataset
    return: NONE (draw histogram)
    """
    path_filout = pathDirectory.result(name_dataset + "/dataset" ) + in_keys.replace (" ","_")
    filout = open (path_filout, "w")
    for PDB_ID in dico_dataset.keys () :
        filout.write (str(dico_dataset[PDB_ID][in_keys]) + "\n")
    filout.close ()
    
    runOtherProg.runRscriptHisto(path_filout, in_keys.replace (" ","_"))
    


def LDAmodel(dico_descriptors_train, dico_descriptors_test, dico_test_correcting, list_descriptor_ttest_signif, path_file_global_train, path_file_global_test, value_cor, path_directory, path_descriptor_signif) :
    
    path_result = pathDirectory.generatePath(path_directory + "validation/")
    
    path_file_train = writeFiles.specificDescriptor(dico_descriptors_train, list_descriptor_ttest_signif, path_result + "validation_train" )
    path_file_test = writeFiles.specificDescriptor(dico_descriptors_test, list_descriptor_ttest_signif, path_result + "validation_test" )
    path_file_test_corected = writeFiles.globalDescriptors(dico_test_correcting, path_result + "validation_test_corrected" )
    
    runOtherProg.ldaPrediction (path_file_train, path_file_test, path_file_global_train, path_file_global_test, path_file_test_corected, value_cor, path_result, path_descriptor_signif)
    
    return pathDirectory.searchModel (path_result)
    
def SVM (dico_descriptors, list_descriptors, path_directory, file_out,  path_file_global_descriptor, name_dataset, train_test_option = 1) :
    """
    SVM predict druggable and non druggable
    args: -> dictionary with descriptor
          -> list descriptors
          -> path directory
          -> name filout
          -> path file global descriptor
          -> name dataset
        option train/test run and elimcor
    return: -> NONE
    """
    
    path_filout = path_directory + file_out

    # train and test
    if train_test_option == 1 : 
        
        path_files_train_test = writeFiles.specificDescriptorbyData (dico_descriptors, list_descriptors, path_filout)
        
        ###############
        #   RUN LDA   #
        ###############
        # run SVM train, test and LOO
        runOtherProg.SVM (path_filin1 =  path_files_train_test [0], path_filin2 =  path_files_train_test [1], path_filout = path_filout)
    else :
        
        ##############
        # Write file #
        ##############
        if list_descriptors == "global" :
            path_file_every_pocket = writeFiles.globalDescriptors(dico_descriptors, path_filout + "_loo")
        else : 
            path_file_every_pocket = writeFiles.specificDescriptor(dico_descriptors, list_descriptors, path_filout + "_loo" )
        
        ###########
        # Run SVM #
        ###########
        runOtherProg.SVM (path_filin1 =  path_file_every_pocket, path_filin2 =  "0", path_filout = path_filout)
   

def scoreBadPredict (list_predicted, dico_dataset_train, dico_dataset_test) : 
    """
    Draw plot with druggability probability
    args: -> list descriptor
          -> dico dataset train
          -> dico dataset test
    return: NONE
    """
    
    path_file_global = list_predicted[0] + "_global"
    filout = open (path_file_global, "w")
    
    
    for file_predicted in list_predicted : 
        if search ("train", file_predicted) :
            path_file_protein = dataSet.predictedByTypePocket (file_predicted, dico_dataset_train) 
            runOtherProg.histNameProtein (path_file_protein)
            filin = open (path_file_protein, "r")
            filin_read = filin.read ()
            filout.write (filin_read)
            filin.close ()
        elif search ("test", file_predicted) : 
            path_file_protein = dataSet.predictedByTypePocket (file_predicted, dico_dataset_test) 
            runOtherProg.histNameProtein (path_file_protein)
            filin = open (path_file_protein, "r")
            filin_read = filin.read ()
            filout.write (filin_read)
            filin.close ()
    
    filout.close ()
    runOtherProg.histNameProtein (path_file_global)
    
    
      
def correlationBetween2Descriptor(dico_descriptors, list_descriptor ) : 
    """
    correlation between two descriptors
    args: -> dictionary descriptors
          -> name descriptors 1
          -> name descriptors 2
    return: NONE draw interest plot
    """
    
    dir_result = pathDirectory.result("correlation")
    path_filout_desc = writeFiles.specificDescriptor(dico_descriptors, list_descriptor, dir_result + list_descriptor[0] +list_descriptor[1] )
    runOtherProg.correlationDescriptor (path_filout_desc)
    

def ACPTwoDatasetDescriptor (dico_descriptor1, dico_descriptor2, list_descriptor, path_dir_result, correspondance_file = "0") : 
    
    path_file_data1 = path_dir_result + "1"
    path_file_data2 = path_dir_result + "2"
    path_file_result = path_dir_result + "ACP_2dataset"
    print dico_descriptor1
    path_file_color = writeFiles.colorACPFile (dico_descriptor1)
    
    if list_descriptor == "global" : 
        writeFiles.globalDescriptors(dico_descriptor1, path_file_data1, debug = 0)
        writeFiles.globalDescriptors(dico_descriptor2, path_file_data2)
    else : 
        writeFiles.specificDescriptor(dico_descriptor1, list_descriptor, path_file_data1, debug = 0 )
        writeFiles.specificDescriptor(dico_descriptor2, list_descriptor, path_file_data2, debug = 0 )
    
    runOtherProg.ACPDataset (path_file_data1, path_file_data2, path_file_result, correspondance_file)

    
    

def histogramFonctionRMSD (dico_descriptor_apo, dico_descriptor_holo, path_file_RMSD, path_dir_result, descriptor) : 
    
    # in here because file write in main. It is not same convention
    
    path_file_descriptor_apo = path_dir_result + "temp_histapo"
    path_file_descriptor_holo = path_dir_result + "temp_histholo"

    writeFiles.specificDescriptor(dico_descriptor_apo, descriptor, path_file_descriptor_apo, debug = 0 )
    writeFiles.specificDescriptor(dico_descriptor_holo, descriptor, path_file_descriptor_holo, debug = 0 )
    
    
    runOtherProg.histoApoHoloDescRMSD(path_file_descriptor_apo, path_file_descriptor_holo, path_file_RMSD, path_dir_result, descriptor)
    
    
    
    
    

def RMSDPockets (dico_dataset, path_dir_result, name_dataset) : 
    
    path_filout = path_dir_result + "RMSDPocket"
    filout = open (path_filout, "w")
    
    for PDB_ID in dico_dataset.keys () : 
        # apo
        if dico_dataset[PDB_ID]["Type structure"] == "apo structure" : 
            
            # apo
            path_dir_descriptor = pathDirectory.descriptor(dir_in = name_dataset + "/Fpocket/" + PDB_ID)
            # directory of PDB_ID
            list_dir = os.listdir(path_dir_descriptor)
#            flag = 0
            for name_dir in list_dir : 
                if search("pocket", name_dir) :
#                    flag = 1 
                    path_pocket_apo_res = pathDirectory.searchPocketResFile(path_dir_descriptor + name_dir + "/")
                    path_pocket_apo_atom = pathDirectory.searchPocketFile(path_dir_descriptor + name_dir + "/")
                    nb_atom_pocket_apo = len(parsePDB.loadCoordSectionPDB(path_pocket_apo_atom[0]))
                    pocket_apo_res = PDB(path_pocket_apo_res)
                    nb_apo_res = len (pocket_apo_res.aaseq())
                    break
#            if flag == 0 : continue
            
            # holo
            PDB_ID_holo = dico_dataset[PDB_ID]["PDB holo"][0]
            print PDB_ID_holo, "*******"
#            print PDB_ID, PDB_ID_holo
            path_dir_descriptor_holo = pathDirectory.descriptor(dir_in = name_dataset + "/Fpocket/" + PDB_ID_holo)
            list_dir_holo = os.listdir(path_dir_descriptor_holo)
            print list_dir_holo
#            flag = 0
            for name_dir_holo in list_dir_holo : 
                if search("pocket", name_dir_holo) : 
#                    flag = 1
                    path_pocket_holo_res = pathDirectory.searchPocketResFile(path_dir_descriptor_holo+ name_dir_holo + "/")
                    print path_pocket_holo_res
                    path_pocket_holo_atom = pathDirectory.searchPocketFile(path_dir_descriptor_holo + name_dir_holo + "/")
                    nb_atom_pocket_holo = len(parsePDB.loadCoordSectionPDB(path_pocket_holo_atom[0]))
                    pocket_holo_res = PDB(path_pocket_holo_res)
                    nb_holo_res = len (pocket_holo_res.aaseq())
                    break
#            if flag == 0 : continue
            
            
            print path_pocket_holo_atom[0], path_pocket_apo_atom[0]
            # superimposition
            list_out_file = runOtherProg.runTMalign(path_pocket_apo_res, path_pocket_holo_res, path_dir_result)
            RMSD = superposeStructure.retrieveRMSDFileTMalign(list_out_file[-1])
            filout.write (str(PDB_ID) + "\t" + str(PDB_ID_holo) + "\t" + str(RMSD) + "\t" + str (nb_apo_res) + "\t" + str (nb_holo_res) + "\t" + str (nb_atom_pocket_apo) + "\t" + str (nb_atom_pocket_holo) + "\t" +pocket_apo_res.aaseq() + "\t" +pocket_holo_res.aaseq() + "\n")
            
    filout.close ()
    
    

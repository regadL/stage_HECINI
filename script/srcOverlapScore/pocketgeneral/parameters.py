"""
BORREL Alexandre
04-2012
run analysis
"""
# personal modules
import pathDirectory
import parseOutLDA
import runOtherProg
import tool
import writeFiles

# global modules
from numpy import arange
import os
from re import search, sub
from itertools import combinations
from copy import deepcopy



def optimizeSVM (path_dir, path_file_global, path_file_train, path_file_test, begin = 0.4, end = 1.01, step = 0.01):
    
    path_dir_result = pathDirectory.generatePath(path_dir + "ElimCor/")
    
    path_filout_sp_sen = path_dir_result + "best_valueSpSen"
    path_filout_acc = path_dir_result + "best_acc"
    filout_sensibility_specificite = open (path_filout_sp_sen, "w")
    filout_acc = open (path_filout_acc, "w")
    for elimcor in arange (begin, end, step) :
        print elimcor, "-> test elimcor"
        path_file_result = runOtherProg.SVM(path_file_global, path_file_train, path_file_test, elimcor=elimcor, debug=1)
        
        try:
            nb_descriptor, acc_loo, acc_train, acc_test = parseOutLDA.retrieveAllAccuracy(path_file_result)
            filout_acc.write (str (elimcor) + "\t" + str (nb_descriptor) + "\t" + str (acc_loo) + "\t" + str (acc_train) + "\t" + str(acc_test) + "\n")
        except:
            pass
        
        try :
            nb_descriptor, accuracy, precision, recall, sensibility, specificity = parseOutLDA.retrieveQualityLOO(path_file_result)
            filout_sensibility_specificite.write (str (elimcor) + "\t" + str(nb_descriptor) + "\t" + str(accuracy) + "\t" + str(precision) + "\t" + str (recall) + "\t" + str (sensibility) + "\t" + str (specificity) + "\n")
        except :
            pass
        
        os.system ("rm " + path_file_result)
        
    filout_sensibility_specificite.close ()
    filout_acc.close ()
    runOtherProg.RplotElimcor(path_filout_sp_sen, 1)
    runOtherProg.RplotElimcor(path_filout_acc, 0)
    
    return bestElimAccuracy (path_filout_acc)




def optimizeLDA (path_dir, path_descriptors, begin = 0.3, end = 1.01, step = 0.01):
    
    path_dir_result = pathDirectory.generatePath(path_dir + "ElimCor/")
    
    path_filout_sp_sen = path_dir_result + "best_valueSpSen"
    path_filout_acc = path_dir_result + "best_acc"
    filout_sensibility_specificite = open (path_filout_sp_sen, "w")
    filout_acc = open (path_filout_acc, "w")
    for elimcor in arange (begin, end, step) :
        runOtherProg.lda(path_descriptors + "_loo", path_filin2=0, path_filout=path_descriptors + "_" + str (elimcor), leave_one_out=1, name_barplot="0", elimcor=elimcor, graph=0, debug=1)
        runOtherProg.lda(path_descriptors + "_train", path_filin2=path_descriptors + "_test", path_filout=path_descriptors + "_" + str (elimcor), leave_one_out=0, name_barplot="0", elimcor=elimcor, graph=0, debug=1)
        
        try:
            nb_descriptor, acc_loo, acc_train, acc_validation = parseOutLDA.retrieveAllAccuracy(path_descriptors + "_" + str (elimcor))
            filout_acc.write (str (elimcor) + "\t")
            filout_acc.write (str (nb_descriptor) + "\t" + str (acc_loo) + "\t" + str (acc_train) + "\t" + str(acc_validation) + "\n")
        except:
            pass
        
        try :
            nb_descriptor, accuracy, precision, recall, sensibility, specificity = parseOutLDA.retrieveQualityLOO(path_descriptors + "_" + str (elimcor))
            filout_sensibility_specificite.write (str (elimcor) + "\t")
            filout_sensibility_specificite.write (str(nb_descriptor) + "\t" + str(accuracy) + "\t" + str(precision) + "\t" + str (recall) + "\t" + str (sensibility) + "\t" + str (specificity) + "\n")
        except :
            pass
        
        
        os.system ("rm " + path_descriptors + "_" + str (elimcor))
        
    filout_sensibility_specificite.close ()
    filout_acc.close ()
    runOtherProg.RplotElimcor(path_filout_sp_sen, 1)
    runOtherProg.RplotElimcor(path_filout_acc, 0)
    
    return bestElimAccuracy (path_filout_acc)


def bestElimCorSpeSen (path_filin):
    
    filin = open (path_filin, "r")
    list_lines = filin.readlines ()
    filin.close ()
    
    min_value = 100
    out_value = 1
    for line_file in list_lines : 
        element_line = line_file.split ("\t")
        elimcor = float (element_line[0])
        specificity = float (element_line[-1])
        sensibility = float (element_line[-2])
        ecart = (1-specificity) + (1-sensibility)
        if ecart < min_value :
            out_value = elimcor
            min_value = ecart
    
    return out_value
    
def bestElimAccuracy (path_filin):
    
    filin = open (path_filin, "r")
    list_lines = filin.readlines ()
    filin.close ()
    
    min_value = 100
    out_value = 1
    for line_file in list_lines : 
        element_line = line_file.split ("\t")
        elimcor = float (element_line[0])
        acc_loo = float (element_line[3])
        acc_train = float (element_line[3])
        acc_test = float (element_line[4])
        
#        if acc_loo == 1 : 
#            acc_loo = 0 
#        if acc_test == 1 : 
#            acc_test = 0
#        if acc_train == 1 : 
#            acc_train = 0
        
        ecart = (1-acc_train) * 0.8 + (1-acc_test)* 0.7  + (1-acc_loo)*0.5
        if ecart < min_value :
            out_value = elimcor
            min_value = ecart
    
    return out_value

def onlyACCLoo (path_filin):
    """
    Select best elimcor for up accuracy
    args: -> path filin
    return: -> best elimcor
    """
    
    acc_temp = 0
        
    for elm in  arange (0.3, 1.1, 0.1) : 
        runOtherProg.lda(path_filin, path_filin2=0, path_filout=path_filin + "_loo", leave_one_out=1, name_barplot="0", elimcor=elm, graph=0, debug=0)
        try : nb_descriptor, accuracy, precision, recall, sensibility, specificity= parseOutLDA.retrieveQualityLOO(path_filin + "_loo")
        except : continue
        if accuracy > acc_temp : 
            out_elim = elm
            acc_temp = accuracy
    
    return out_elim
    
    
def retrieveDescriptorForLDAModel (dico_desc, path_file_train, path_file_test, path_dir, type_selection = "AllModel", begin = 0.01, end = 1.01, step = 0.01, nb_descriptor = 0, nb_model = 100, debug = 0):
    """
    Retrieve list of descriptors 
    args : -> path file train
           -> path file test
    return: -> list descriptors
    """
    # lot of fonction -> elimcor (not use now)
    # selection with all model
 
    ###### short cut 
    path_file_descriptor_ACC = path_dir + "tempSelected"
#     path_file_descriptor_ACC = runOtherProg.descriptorSelectionLDA (path_file_train, path_file_test,  path_dir + "tempSelected", nb_descriptor)

    if type_selection == "AllModel" : 
        list_descriptor = selectBestModelAmongBestModel (dico_desc, path_file_descriptor_ACC, nb_model)
 
    return list_descriptor
 
def retrieveDescriptorForSVMModel (dico_descriptor, path_dir, nb_descriptor = 3, debug = 0):
    """
    Retrieve list of descriptors 
    args : -> path file train
           -> path file test
    return: -> list descriptors
    """
    
    l_descriptors = tool.listDescriptor(dico_descriptor["Druggable"])
    print l_descriptors
    
    print len (l_descriptors)
    l_combi_desc = list(combinations(l_descriptors, nb_descriptor))
    print len (l_combi_desc)
    
    i_combi = len (l_combi_desc) - 2
    nb_combi = len (l_combi_desc)
    
    while i_combi < nb_combi : 
        
        l_p_file = writeFiles.specificDescriptorbyData(dico_descriptor, list(l_combi_desc[i_combi]), path_dir + str(i_combi), debug = 0)
        runOtherProg.descriptorSelectionSVM(l_p_file[0], l_p_file[1], str (nb_descriptor) + "_" + str(i_combi) + ".svm")
        
        i_combi = i_combi + 1
        
    
    
    # lot of fonction -> elimcor (not use now)
    # selection with all model
 
    ###### pass line ######
#    path_file_descriptor_ACC = path_dir + "tempSelected"

#     path_file_descriptor_ACC = runOtherProg.descriptorSelectionSVM (path_file_train, path_file_test,  path_dir + "tempSelected", nb_descriptor)

#     if type_selection == "AllModel" : 
#         list_descriptor = selectBestModelAmongBestModel (path_file_descriptor_ACC)
 
#     print list_descriptor
#     return list_descriptor
 
 
#    AFTER with elimcor   
#    filout_acc = open (path_prefix + "_selectDesc", "w")
#    filout_desc = open (path_prefix + "_Descriptor", "w")
    
#    for value_cor in arange (begin, end, step) : 
#        path_file_descriptor_ACC = runOtherProg.descriptorSelection (path_file_train, path_file_test, value_cor, path_prefix + "_tempACC")
#        list_descriptor , acc  = tool.parseSelectDescriptor (path_file_descriptor_ACC)
#        filout_acc.write (str(value_cor) + "\t" + str (len (list_descriptor) - 1) + "\t" + str (acc[2]) + "\t" + str (acc[0]) + "\t" +  str (acc[1]) + "\n")
#        filout_desc.write (str(value_cor) + "\t" + "\t".join(list_descriptor) + "\n")
#    filout_acc.close ()
#    filout_desc.close ()
#    
#    runOtherProg.RplotElimcor(path_prefix + "_selectDesc", 0)
#    best_cor_value = bestElimAccuracy(path_prefix + "_selectDesc")
#    
#    path_file_descriptor_ACC = runOtherProg.descriptorSelection (path_file_train, path_file_test, best_cor_value, path_prefix + "_bestACC")
#    list_descriptor_out , acc  = tool.parseSelectDescriptor (path_file_descriptor_ACC)
#    
#    list_descriptor_out.remove ("drugg")
#    if debug :
#        print ("************")
#        print (best_cor_value)
#        print (len(list_descriptor_out))
#        print ("************")    
    
#    return list_descriptor


def countDescriptorFile (dico_desc, path_filin, number_model_selected) :
    """
    Count descriptor in list of model
    args -> path file with descriptor
         -> number of best model selected
    return: -> dictionary count
            -> combinaison descriptor
            -> best model
    """
    p_filout = os.path.dirname (path_filin) + "/acc_mcc"
    filout = open (p_filout, "w")
    l_d_model = []
    
    criterion_mcc = 0
    filin = open (path_filin, "r")
    filin_read = filin.read ()
    list_model = filin_read.split ("[1] \"**********************************************************************\"\n") 

    for model in list_model : 
        l_desc_acc_mcc_spse = retrieveDescriptorAccModel(model.split ("\n"))
        if l_desc_acc_mcc_spse : 
            d_model = {}
            l_desc = l_desc_acc_mcc_spse [0]
            l_acc = l_desc_acc_mcc_spse [1]
            l_mcc = l_desc_acc_mcc_spse [2]
            l_sp_se = l_desc_acc_mcc_spse[3]
            l_des_implication = l_desc_acc_mcc_spse[4]
            
            if l_mcc[0] == 'NaN' or l_mcc[1] == 'NaN' or l_mcc[2] == 'NaN' : 
                continue
            # implement model
            d_model["acc LOO"] = float (l_acc[0])
            d_model["acc TRAIN"] = float (l_acc[1])
            d_model["acc TEST"] = float (l_acc[2])
            
            d_model["descriptors"] = l_desc
            d_model["implication"] = l_des_implication
            
            d_model["mcc LOO"] = float (l_mcc[0])
            d_model["mcc TRAIN"] = float (l_mcc[1])
            d_model["mcc TEST"] = float (l_mcc[2])
            
            d_model["se LOO"] = float (l_sp_se[0])
            d_model["se TRAIN"] = float (l_sp_se[2])
            d_model["se TEST"] = float (l_sp_se[4])
            
            d_model["sp LOO"] = float (l_sp_se[1])
            d_model["sp TRAIN"] = float (l_sp_se[3])
            d_model["sp TEST"] = float (l_sp_se[5])
            
            # criteria
            d_model["criterion mcc"] = float (l_mcc[0]) + float (l_mcc[2]) # MCC with LOO
            d_model["criterion acc"] = float (l_acc[0]) + float (l_acc[1])
            
            l_d_model.append (d_model)
            
            # control rates 
            filout.write ("\t".join(l_acc) + "\t" + "\t".join(l_mcc) + "\n")
            
            # select best model, win one loop
            if d_model["criterion mcc"] > criterion_mcc : 
                best_model = l_desc
                #tool.percentageR(l_desc)
                criterion_mcc = d_model["criterion mcc"]

    
    # histogram to control the score selected -> distribution
    filout.close ()
    runOtherProg.plotQualityMCCACC (p_filout)

    # object programation -> change directly the list of model -> limit the number
    selectBestModel (dico_desc, l_d_model, os.path.dirname (path_filin), number_model_selected)
    
    # dico count            
    dico_count = countDescriptor (l_d_model)
    dico_count_combi = countDescriptorCombi (l_d_model)
    dico_desc_ranking = countRakingDesc (l_d_model)    
#     tool.printDico(dico_count)
#     print (dico_count_combi)
    return [dico_count, dico_count_combi, dico_desc_ranking, best_model]


def retrieveDescriptorAccModel (list_lines) : 

    nb_lines = len (list_lines)

    i=0
    while i< nb_lines : 
        if search ("\[1\] \"descriptor\"", list_lines[i]) :
            list_descriptor = tool.formatLinesDescriptor (list_lines[i + 1])
            if not search ("^\[1\]",list_lines[i+2]) : 
                list_descriptor = list_descriptor + tool.formatLinesDescriptor(list_lines[i+2])
                l_des_implication = sub("[ ]{2,}", " ",list_lines[i+3].strip())
                l_des_implication = l_des_implication.split (" ")[1:]
                acc_global = tool.formatLinesAcc(list_lines[i+5])
            else :
                acc_global = tool.formatLinesAcc(list_lines[i+4])
                l_des_implication = sub("[ ]{2,}", " ",list_lines[i+2].strip())
                l_des_implication = l_des_implication.split (" ")[1:]
                
                
            mcc_global = tool.formatLinesAcc(list_lines[-2])
            spse_global = tool.formatLinesAcc(list_lines[-4])
            
            return [list_descriptor, acc_global, mcc_global, spse_global, l_des_implication]
        i = i + 1
    


def countDescriptor (l_model) : 

    d_count = {}
    
    for model in l_model : 
        for desc in model["descriptors"] : 
            if not desc in d_count.keys () : 
                d_count[desc] = 1
            else : 
                d_count[desc] =  d_count[desc] + 1
    


    return d_count


def countDescriptorCombi (l_model) : 
    
    dico_count = {}
    
    for model in l_model : 
        list_descriptor = model["descriptors"]
        nb_descriptor = len (list_descriptor)
        
        i = 0
        while i < nb_descriptor - 1 : 
            #print i, "i"
            j = i + 1
            while j< nb_descriptor : 
                #print j, "j"
                if not list_descriptor[i] in dico_count.keys () :
                    if not list_descriptor[j] in  dico_count.keys () :
                        dico_count[list_descriptor[i]] = {}
                        dico_count[list_descriptor[i]][list_descriptor[j]] = 1
                    else :
                        if not list_descriptor[i] in dico_count[list_descriptor[j]].keys () : 
                            dico_count[list_descriptor[j]][list_descriptor[i]] = 1    
                        else :             
                            dico_count[list_descriptor[j]][list_descriptor[i]] = dico_count[list_descriptor[j]][list_descriptor[i]] + 1
                else : 
                    if not list_descriptor[j] in dico_count[list_descriptor[i]].keys () :
                        dico_count[list_descriptor[i]][list_descriptor[j]] = 1
                    else :                
                        dico_count[list_descriptor[i]][list_descriptor[j]] = dico_count[list_descriptor[i]][list_descriptor[j]] + 1
                j = j + 1
            i = i + 1

    return dico_count
                

def countRakingDesc (l_d_model)  : 
    
    d_out = {}
    d_out["descriptors"] = []
    nb_desc = len (l_d_model[0]["descriptors"])
    
    print nb_desc, "nb_des"
    print l_d_model[0]["descriptors"]
    
    i = 0
    while i < nb_desc : 
        d_out[i] = {}
        i = i + 1
    
    for model in l_d_model : 
        l_sort = deepcopy(model["implication"])
        l_sort.sort (reverse=True)
        
        i = 0
        while i < nb_desc : 
            print i, l_sort
            i_val = model["implication"].index (l_sort[i])
                
            if not model["descriptors"][i_val] in d_out[i].keys () : 
                d_out[i][model["descriptors"][i_val]] = 1
            else : 
                d_out[i][model["descriptors"][i_val]] = d_out[i][model["descriptors"][i_val]] + 1
 
            if not model["descriptors"][i_val] in d_out["descriptors"] : 
                d_out["descriptors"].append (model["descriptors"][i_val])
            else : 
                pass
 
            i = i + 1
            
             
    return d_out
    
    
    



def selectBestModelAmongBestModel (dico_desc, path_file_every_best_model, nb_best_selected) : 
    
    # open file with model and select first model
    table_count =  countDescriptorFile (dico_desc, path_file_every_best_model, nb_best_selected)
    
    path_file_count = writeFiles.writeCount(table_count[0], path_file_every_best_model + str (nb_best_selected) + "_countTable")
    path_file_count_combi = writeFiles.writeCountCombi(table_count[1], path_file_every_best_model + str (nb_best_selected) + "_countCombiTable")
    
    p_file_implication = writeFiles.descImplicationLDA (table_count[2], path_file_every_best_model + str (nb_best_selected) + "_countImplication")
    
    runOtherProg.histNameProtein(path_file_count)
    runOtherProg.histNameProtein(path_file_count_combi)
    runOtherProg.pieDescriptor (path_file_count)
    runOtherProg.pieDescriptor (p_file_implication)
    
    # best model in output table count 
    return table_count[3]
    
    
    
def  selectBestModel (dico_descriptor, l_model, p_dir, number_model_selected, debug = 0) : 
    
    
    # criterion selected
    l_criterion = []
    # histogram score MCC
    p_filout_histo = p_dir + "/score_MCC_hist_" + str (number_model_selected)
    p_filout_model = p_dir + "/best_" + str(number_model_selected) + "_model"
    p_filout_score = p_dir + "/score_criterion_" + str(number_model_selected)
    p_dir_model = pathDirectory.generatePath(p_dir + "/BestModels" + str (number_model_selected)+"/")
    
    filout_histo = open(p_filout_histo, "w")
    filout_model = open(p_filout_model, "w")
    filout_score = open (p_filout_score, "w")
    
    nb_model = len (l_model)
    if debug : print len (l_model)
    i = 0
    while i < nb_model : 
        l_criterion.append( l_model[i]["criterion mcc"])
        i = i + 1
    
    if debug : print len(l_model)
    
    l_criterion.sort(reverse=True)
    print l_criterion[0:13]
    

    if len (l_criterion) > number_model_selected : 
        criterion = l_criterion[number_model_selected+1]
    else : 
        criterion = l_criterion[-1]
    
    i = 0
    nb_model = len (l_model)
    
    if debug :print criterion, "criterion"
    if debug : print len (l_model) 
    
    while i < nb_model : 
        if l_model[i]["criterion mcc"] < criterion : 
            del l_model[i]
            nb_model = nb_model - 1
            continue
        else :
            # write file with n best models
            filout_model.write ("************************************************\n")
            filout_model.write ("----".join(l_model[i]["descriptors"]) + "\n")
            filout_model.write ("************--ACCURACY--**********\n")
            filout_model.write (str(l_model[i]["acc LOO"]) + "---" + str(l_model[i]["acc TRAIN"]) + "---" + str(l_model[i]["acc TEST"]) + "\n")
            filout_model.write ("************--SP SE--**********\n")
            filout_model.write (str(l_model[i]["sp LOO"]) + "---" + str(l_model[i]["se LOO"]) + "---" + str(l_model[i]["sp TRAIN"]) + "---" + str(l_model[i]["se TRAIN"]) + "---" + str(l_model[i]["sp TEST"]) + "---" + str(l_model[i]["se TEST"]) + "\n")
            filout_model.write ("************--MCC--**********\n")
            filout_model.write (str(l_model[i]["mcc LOO"]) + "---" + str(l_model[i]["mcc TRAIN"]) + "---" + str(l_model[i]["mcc TEST"]) + "\n")
            
            filout_histo.write (str (l_model[i]["criterion mcc"]) + "\t" + str(l_model[i]["mcc TEST"]) + "\n")
            filout_score.write ("M" + str (i + 1) + "\t" + str (l_model[i]["criterion mcc"]) + "\t" + str(l_model[i]["mcc LOO"]) + "\t" + str(l_model[i]["mcc TEST"]) + "\n" )
            # make model in folder
            path_desc = writeFiles.specificDescriptorbyData(dico_descriptor, l_model[i]["descriptors"], p_dir_model + str(i))[0]
            runOtherProg.runSaveModelLDA(path_desc)
            
            i = i + 1 
    
    if debug : print len (l_model)
    filout_histo.close ()
    filout_model.close ()
    filout_score.close ()
    
    # control score distribution
    runOtherProg.runRscriptHisto (p_filout_histo, "MCC_criterion", brk = 10, order = number_model_selected)
    matrixImageModel (l_model, p_dir + "/image_model_" + str (number_model_selected))
    runOtherProg.matrixCorr (p_dir + "/image_model_" + str (number_model_selected), 2, debug=1)
    runOtherProg.matrixImage (p_dir + "/image_model_" + str (number_model_selected), 2, debug=1)
    runOtherProg.plotScoreCriterion (p_dir + "/score_criterion_" + str(number_model_selected))
    
    
def matrixImageModel (l_model, p_filout):
    l_desc_head = []
    for model in l_model : 
        for desc in model ["descriptors"] : 
            if not desc in l_desc_head : 
                l_desc_head.append (desc)
                
    filout = open (p_filout, "w")
    filout.write ("\t".join (l_desc_head) + "\n")
    
    l_write = []
    i = 1
    for model in l_model :
        l_write = [] 
       
        for desc_head in l_desc_head :
            flag = 0
            for desc_model in model["descriptors"] :  
                if desc_model == desc_head : 
                    l_write.append ("1")
                    flag = 1
            if flag ==0 : 
                l_write.append ("0")
        
        filout.write ("M" + str (i) + "\t" + "\t".join (l_write) + "\n")
        i = i + 1
    
    filout.close ()
            
            
            
    
    
    
           
    


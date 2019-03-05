import dataSet
import writeFiles



def qualityFpocket (pocket_retrieve_type, name_dataset, path_filout):
    
    dico_fpocket_global = dataSet.loadDatasetWithFpocketDrugScore (name_dataset, pocket_retrieve_type)
    print dico_fpocket_global

    
    filout = open (path_filout, "w")
    filout.write ("------Fpocket Druggability Score-------\n")
    filout.write ("------GLOBAL-------\n")
    dico_rate_global = FpocketRetrieveRate(dico_fpocket_global)
    quality(dico_rate_global) # append quality rate
    writeFiles.rateQualityFile (dico_rate_global, filout)
    
    if name_dataset == "krasowski" : 
        filout.write ("------TRAINNING-------\n")
        dico_Fpocket_train = dataSet.reduceDicoDataset (dico_fpocket_global, "data", "t")
        dico_rate_train = FpocketRetrieveRate(dico_Fpocket_train)
        quality(dico_rate_train) # append quality rate
        writeFiles.rateQualityFile (dico_rate_train, filout)
    
        filout.write ("------VALIDATION-------\n")
        dico_Fpocket_test = dataSet.reduceDicoDataset (dico_fpocket_global, "data", "v")
        dico_rate_test = FpocketRetrieveRate(dico_Fpocket_test)
        quality(dico_rate_test) # append quality rate
        writeFiles.rateQualityFile (dico_rate_test, filout)  
    
    filout.close ()



def quality (rate):
    
    rate["acc"] = ((rate["TP"] + rate["TN"])/(rate["TP"] + rate["TN"] + rate["FN"] + rate["FP"]))
    rate["pr"] = ((rate["TP"])/(rate["TP"] + rate["FP"]))
    rate["se"] = ((rate["TP"])/(rate["TP"] + rate["FN"]))
    rate["sp"] = ((rate["TN"])/(rate["TN"] + rate["FP"]))
    
    


def FpocketRetrieveRate (dico_fpocket, debug = 1):
   
    dico_out = {}
    dico_out["TP"] = 0.0
    dico_out["FP"] = 0.0
    dico_out["TN"] = 0.0
    dico_out["FN"] = 0.0
    
    
    for PDB in dico_fpocket.keys () : 
        print PDB
        print dico_fpocket[PDB]
        if dico_fpocket[PDB]["druggability"] == "n" :
            
            if not "Drug Score" in dico_fpocket[PDB].keys ()  : 
                continue
             
            if float(dico_fpocket[PDB]["Drug Score"]) <= 0.5 :
                dico_out["TN"] = dico_out["TN"] + 1
            else : 
                dico_out["FP"] = dico_out["FP"] + 1
        else : 
            if  not "Drug Score" in dico_fpocket[PDB].keys():
                continue 
            
            if float(dico_fpocket[PDB]["Drug Score"]) <= 0.5 :
                dico_out["FN"] = dico_out["FN"] + 1
            else : 
                dico_out["TP"] = dico_out["TP"] + 1 
                
    
    if debug : 
        print dico_out
    
    return dico_out






        
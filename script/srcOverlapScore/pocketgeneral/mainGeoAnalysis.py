from os import path, listdir
from re import search
import dataSet
import pathDirectory
import checkPocket
import loadDescriptors
import descriptor
import tool
import writeFiles
import runOtherProg
import overlapPocket



path_file_krasowski = pathDirectory.dataSet() + "krasowski_2011_dataset.txt"




def descriptorCompute (estimator_type, name_dataset):
    
    # dataset
    if name_dataset == "krasowski" : 
        l_pdb = dataSet.formatFileDataSetKrasowki2011(path_file_krasowski).keys ()
    else : 
        pass
    
    p_dir_dataset = pathDirectory.dataSet(name_dataset)
    
    
    if estimator_type == "DOGSITE" : 
        l_dir_pocket = checkPocket.retrievePocketDogSite (p_dir_dataset, l_pdb, debug = 1) # move result
        # run descriptor geo -> radi
        
        print l_dir_pocket
        
        for dir_pocket in l_dir_pocket : 
            l_p_file = listdir(dir_pocket)
            
            for p_file in l_p_file : 
                if search(".pdb", p_file) : 
                    path_pocket_atom = dir_pocket + p_file
                    ##########################################
                    # Check PDB end file for RADI, END final #
                    ##########################################
#                     tool.checkENDFinalLinePDBfile(path_pocket_atom)
#                     tool.checkHEADERinitialLinePDBfile(path_pocket_atom)
#                     descriptor.radi (dir_pocket + p_file, dir_pocket + p_file[:-4] , "none", "dogsite")

    elif estimator_type == "Fpocket" : 
#         for PDB in l_pdb : 
#             p_dir_pocket = p_dir_dataset + PDB + "/protein_out/pockets/"
#             
#             l_p_file = listdir(p_dir_pocket)
#             
#             for p_file in l_p_file : 
#                 if search("atm.pdb", p_file) : 
#                     path_pocket_atom = p_dir_pocket + p_file
                    ##########################################
                    # Check PDB end file for RADI, END final #
                    ##########################################
#                     tool.checkENDFinalLinePDBfile(path_pocket_atom)
#                     tool.checkHEADERinitialLinePDBfile(path_pocket_atom)
#                     descriptor.radi (p_dir_pocket + p_file, p_dir_pocket + p_file[:-4] , "none", "Fpocket")
            
                
        dico_est = loadDescriptors.FpocketGlobal (p_dir_dataset, l_pdb)
    
    

def histVolume (estimator_type, name_dataset ):

    p_dir_result = pathDirectory.result("GeoStudy")

    if name_dataset == "krasowski" : 
        l_pdb = dataSet.formatFileDataSetKrasowki2011(path_file_krasowski).keys ()
    
    p_dir_dataset = pathDirectory.dataSet(name_dataset)

    if estimator_type == "DOGSITE" : 
        dico_est = loadDescriptors.DogSiteGlobal (p_dir_dataset, l_pdb)
        # -----
        l_desc = dico_est[dico_est.keys ()[0]]
        l_desc_est =  l_desc[l_desc.keys()[0]].keys ()
        # -------
        for desc in l_desc_est :
            try : 
                name_file = desc.replace (" ", "_")
            except :
                name_file = desc
            
            try : 
                name_file = desc.replace ("/", "_")
            except :
                name_file = desc
              
            p_vol_txt = writeFiles.listValue(dico_est, desc, p_dir_result +  estimator_type + name_file + "_DogSite.txt")
            runOtherProg.runRscriptHisto (p_vol_txt, name_file)
        
     
         
    elif estimator_type == "Fpocket" : 
        dico_est = loadDescriptors.FpocketGlobal (p_dir_dataset, l_pdb)
        # -----
        l_desc = dico_est[dico_est.keys ()[0]]
        l_desc_est =  l_desc[l_desc.keys()[0]].keys ()
        # -------
        for desc in l_desc_est : 
            try : 
                name_file = desc.replace (" ", "_")
            except :
                name_file = desc
                
            p_vol_txt = writeFiles.listValue(dico_est,desc, p_dir_result +  estimator_type + name_file + "_Fpocket.txt")
            runOtherProg.runRscriptHisto (p_vol_txt, name_file)
        
        
            
#         print tool.printDico(dico_fpocket)
        
        
                
    dico_radi = loadDescriptors.radiglobal(p_dir_dataset, l_pdb, estimator_type)
    # -------
    l_desc_radi = dico_radi[dico_radi.keys ()[0]]
    l_desc_radi =  l_desc_radi[l_desc_radi.keys()[0]].keys ()
    # -------

    for desc_radi in l_desc_radi : 
        try : 
            name_file = desc_radi.replace (" ", "_")
        except :
            name_file = desc_radi
            
        # radi    
        p_radi_txt = writeFiles.listValue(dico_radi, desc_radi, p_dir_result +  estimator_type + name_file + "_radi.txt")
        runOtherProg.runRscriptHisto (p_radi_txt, name_file)
     




def correlationVolume (name_dataset, type_desc_vol = "", thresholdSO = 0.5):
    """
    Correlation DogSite-Fpocket
    """
    
    p_dir_dataset = pathDirectory.dataSet(name_dataset)
    p_d_result = pathDirectory.result("GeoStudy")
    
    if name_dataset == "krasowski" : 
        l_pdb = dataSet.formatFileDataSetKrasowki2011(path_file_krasowski).keys ()
    
    # load descriptor 
    if type_desc_vol == "radi" : 
        dico_vol_Fpock = loadDescriptors.radiglobal(p_dir_dataset, l_pdb, "Fpocket")
        dico_vol_DogSite = loadDescriptors.radiglobal(p_dir_dataset, l_pdb, "DOGSITE")
        
    else : 
        dico_vol_Fpock = loadDescriptors.FpocketGlobal(p_dir_dataset, l_pdb)
        dico_vol_DogSite = loadDescriptors.DogSiteGlobal(p_dir_dataset, l_pdb)
        
    
    p_filout = p_d_result + "histSOVol_" + type_desc_vol + "_" + str (thresholdSO) + ".txt"
    filout = open (p_filout, "w")
    filout.write ("SO\tFpocket\tDoGSite\n")
    
    for PDB in l_pdb :
        
        p_dir_Fpocket = p_dir_dataset + PDB + "/protein_out/pockets/"
        p_dir_dogsite = p_dir_dataset + PDB + "/DOGSITE/"
        
        l_p_pockfpocket = listdir(p_dir_Fpocket)
        l_p_pockDogsite = listdir(p_dir_dogsite)
        
        for p_pockfpocket in l_p_pockfpocket : 
            if search("atm.pdb", p_pockfpocket) : 
                for p_pockDogSite in l_p_pockDogsite : 
                    if search(".pdb", p_pockDogSite) : 
                        SO = overlapPocket.scoreOverlap(p_dir_dogsite + p_pockDogSite, p_dir_Fpocket + p_pockfpocket)
                        if SO > 0 : 
                            if type_desc_vol != "radi" : 
                                try : filout.write (str(SO) + "\t" + str (dico_vol_Fpock[PDB][p_dir_Fpocket + p_pockfpocket]["Real volume (approximation)"]) + "\t" + str (dico_vol_DogSite[PDB][p_dir_dogsite + p_pockDogSite]["Volume[A^3]"]) + "\n" )
                                except : pass
                            else : 
                                filout.write (str(SO) + "\t" + str (dico_vol_Fpock[PDB][p_dir_Fpocket + p_pockfpocket]["VOLUME_HULL"]) + "\t" + str (dico_vol_DogSite[PDB][p_dir_dogsite + p_pockDogSite]["VOLUME_HULL"]) + "\n" )
#                                 filout.write (str(SO) + "\t" + str (dico_vol_Fpock[PDB][p_dir_Fpocket + p_pockfpocket]["PCI"]) + "\t" + str (dico_vol_DogSite[PDB][p_dir_dogsite + p_pockDogSite]["PCI"]) + "\n" )
    filout.close ()
    
    
    # plot correlation
    runOtherProg.correlationPlotSO (p_filout, thresholdSO)
    
                           
         
def correlationCardRadi (name_dataset, estimator_type):
    
    p_dir_dataset = pathDirectory.dataSet(name_dataset)
    p_d_result = pathDirectory.result("GeoStudy")
    
    if name_dataset == "krasowski" : 
        l_pdb = dataSet.formatFileDataSetKrasowki2011(path_file_krasowski).keys ()
    
    # load descriptor 
    dico_radi = loadDescriptors.radiglobal(p_dir_dataset, l_pdb, estimator_type)
    
    if estimator_type == "Fpocket" : 
        dico_est = loadDescriptors.FpocketGlobal(p_dir_dataset, l_pdb)
#         tool.printDico(dico_est)
    else :
        dico_est = loadDescriptors.DogSiteGlobal(p_dir_dataset, l_pdb)
    
    
    l_desc = dico_radi[dico_radi.keys ()[0]]
    l_desc_radi =  l_desc[l_desc.keys()[0]].keys ()
    
    
    l_desc = dico_est[dico_est.keys ()[0]]
    l_desc_est =  l_desc[l_desc.keys()[0]].keys ()
    
    try : l_desc_est.remove ("Name")
    except : pass
    
     
    p_filout = p_d_result + "matrixDesc_" + estimator_type +".txt"
    filout = open (p_filout, "w")
    
    # header 
    filout.write ("\t".join(l_desc_est) + "\t")
    filout.write ("\t".join(l_desc_radi) + "\n")
    
    i=0
     
    for PDB in l_pdb :
        # DogSite
        if estimator_type == "DOGSITE" : 
            p_dir = p_dir_dataset + PDB + "/DOGSITE/"
            l_p_pock = listdir(p_dir)
            for p_pock in l_p_pock : 
                if search(".pdb", p_pock) : 
                    p_file = p_dir + p_pock
            
                    # try input dictionnary
                    try : 
                        test = dico_est[PDB][p_file][l_desc_est[0]]
                    except :
                        continue
                    
                                        
                    l_val = []    
                    for des_est in l_desc_est : 
                        try : l_val.append (str(dico_est[PDB][p_file][des_est]))
                        except : pass
                    for des_radi in l_desc_radi : 
                        l_val.append (str(dico_radi[PDB][p_file][des_radi]))
                        
                    filout.write ("\t".join(l_val) + "\n")
                    
                    
            
        #Fpocket    
        elif estimator_type == "Fpocket" : 
            p_dir = p_dir_dataset + PDB + "/protein_out/pockets/"
            l_p_pock = listdir(p_dir)
            for p_pock in l_p_pock : 
                if search("atm.pdb", p_pock) : 
                    p_file = p_dir + p_pock
                
                    l_val = []    
                    for des_est in l_desc_est : 
                        l_val.append (str(dico_est[PDB][p_file][des_est]))
                    for des_radi in l_desc_radi : 
                        l_val.append (str(dico_radi[PDB][p_file][des_radi]))
                    filout.write ("\t".join(l_val) + "\n")
        
    filout.close ()
    
    runOtherProg.matrixCorr(p_filout, 6) 
     
     
    
def PCADescriptorGeo (name_dataset, estimator) : 
    
    p_dir_result = pathDirectory.result("GeoStudy")
    p_dir_dataset = pathDirectory.dataSet(name_dataset)
    
    
    if name_dataset == "krasowski" : 
        l_pdb = dataSet.formatFileDataSetKrasowki2011(path_file_krasowski).keys ()
    
    
    if estimator == "Fpocket" : 
        dico_est = loadDescriptors.FpocketGlobal (p_dir_dataset, l_pdb)
        l_desc = ["Number of V. Vertices", "Mean B-factor", "Mean alpha-sphere radius", "Number of apolar alpha sphere", "Real volume (approximation)"]
    elif estimator == "DOGSITE" : 
        dico_est = loadDescriptors.DogSiteGlobal(p_dir_dataset, l_pdb)
        l_desc = ["Ellips b/a", "Ellips c/a", "Depth[A]", "Surface[A^2]", "enclosure", "Volume[A^3]", "Lipo_surf[A^2]"]
#         l_desc = ["Volume[A^3]"]
   
   
    # load radi
    dico_radi = loadDescriptors.radiglobal(p_dir_dataset, l_pdb, estimator)
    l_desc_radi = dico_radi[dico_radi.keys ()[0]]
    l_desc_radi =  l_desc_radi[l_desc_radi.keys()[0]].keys ()
    
#     l_desc_radi = ["PCI", "RADIUS_CYLINDER"]
     
        
    p_filout = p_dir_result + "PCA_" + estimator +".txt"
    filout = open (p_filout, "w")    
    
    # header 
    filout.write ("\t".join(l_desc) + "\t")
    filout.write ("\t".join(l_desc_radi) + "\n")
    
    for PDB in l_pdb : 
        for p_file in dico_est[PDB].keys () : 
            l_val = [] 
             
            # check file exist
            try :   
                t = dico_est[PDB][p_file][l_desc[0]]
                tt = dico_radi[PDB][p_file][l_desc_radi[0]]
            except :
                continue
            
            # descriptor
            for des_est in l_desc : 
                l_val.append (str(dico_est[PDB][p_file][des_est]))
            for des_radi in l_desc_radi : 
                l_val.append (str(dico_radi[PDB][p_file][des_radi]))
            filout.write ("\t".join(l_val) + "\n")
    filout.close ()
    
    #writeFiles.colorACPFile(dico_descriptor) -> modifier la fonction en dure !!!
    runOtherProg.runRACP(p_filout, "ACP_geo_" + estimator, 1)
    
    
    






    
    
          

# histVolume("Fpocket", "krasowski")
# histVolume("DOGSITE", "krasowski")


    
# correlationVolume ("krasowski", "radi")
# correlationVolume ("krasowski")

correlationCardRadi ("krasowski","DOGSITE")
correlationCardRadi ("krasowski","Fpocket")


# PCADescriptorGeo ("krasowski","DOGSITE")
# PCADescriptorGeo ("krasowski","Fpocket")




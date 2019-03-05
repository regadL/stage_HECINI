# personnal module
import StructureStock
import pathDirectory
import writeFiles
import parseDescriptorRadii3
import tool
import descriptor
import parseFpocket

# general module
import os, re

def generalLoadDescriptor (dico_dataset, type_protomol, type_pocket, name_dataset, only_one_pocket = 1, write_file = 1, option_separate = 1, debug = 0):
    """
    Load descriptors for every pockets
    args: -> ditionary with dataset
          -> type protomol (surflexe or Fpocket)
          -> select only PDB with only one pocket (default 1)
          -> write file with descriptors (default 1 -> write)
          -> drugg separation (separate default 1)
    return : NULL write files
    """
    
    dico_descriptor = {}
    dico_descriptor["Druggable"] = {}
    dico_descriptor["Druggable"]["PDB"] = []
    dico_descriptor["Druggable"]["data"] = []
    dico_descriptor["No-Druggable"] = {}
    dico_descriptor["No-Druggable"]["PDB"] = []
    dico_descriptor["No-Druggable"]["data"] = []
    dico_descriptor["data"] = []
    
    list_PDB = dico_dataset.keys ()
    for PDB_ID in list_PDB :
        # pocket gestion directory -> remove file
        list_dir_pocket = pathDirectory.generateListDirPocket( PDB_ID,type_pocket, name_dataset)
        if list_dir_pocket == [] : # bug pocket estimation 
            print PDB_ID, "-> No pocket detected"
            continue 
        
        # append in structure
        if option_separate : 
            if dico_dataset[PDB_ID] ["druggability"] == "d" : 
                dico_load = dico_descriptor["Druggable"]
            else  : 
                dico_load = dico_descriptor["No-Druggable"]
        else : 
            dico_load = dico_descriptor["No-Druggable"]
            
        if "data" in dico_dataset[PDB_ID].keys() :     
            dico_load["data"].append (dico_dataset[PDB_ID]["data"])
        else : pass
        
        # pocket gestion directory
        list_dir_pocket = pathDirectory.generateListDirPocket( PDB_ID,type_pocket, name_dataset)
            
        if list_dir_pocket == [] : # bug pocket estimation 
            print PDB_ID, "-> No pocket detected"
            continue 
        
        if debug :
            print list_dir_pocket,dico_dataset[PDB_ID] ["druggability"] 
        
        if only_one_pocket == 1 and len (list_dir_pocket) > 1 : # check number pocket
            list_dir_pocket = [list_dir_pocket[0]]
            #del dico_dataset[PDB_ID]
            #continue
        for dir_pocket in list_dir_pocket : 
#            loadDescriptorPocket (dir_pocket, dico_load, "area", type_protomol, type_pocket) 
            loadDescriptorPocket (dir_pocket, dico_load, "atomic", type_protomol, type_pocket)
            loadDescriptorPocket (dir_pocket, dico_load, "energy", type_protomol, type_pocket)
            loadDescriptorPocket (dir_pocket, dico_load, "volume", type_protomol, type_pocket)
            loadDescriptorPocket (dir_pocket, dico_load, "radi", type_protomol, type_pocket)
#             loadDescriptorFPocket (dir_pocket, dico_load, type_protomol, type_pocket)
    
        #print dico_load
        dico_load["PDB"].append (PDB_ID)
        tool.appendNAValue (dico_load)
        
        if debug :
            print tool.nbValue(dico_load), "next NA"
            print len (dico_load["PDB"]), dico_load["PDB"]
            print dico_load
    
    #descriptor.ratioCalculation (dico_descriptor) a finir
    if write_file : 
        if option_separate : 
            writeFiles.globalDescriptors(dico_descriptor, pathDirectory.result(name_dataset + "/" + type_pocket) + "global_descriptor_druggable_" + type_protomol +"_" + type_pocket,"Druggable")
            writeFiles.globalDescriptors(dico_descriptor, pathDirectory.result(name_dataset + "/" + type_pocket) + "global_descriptor_nodruggable_" + type_protomol +"_" + type_pocket, "No-Druggable")
        writeFiles.globalDescriptors(dico_descriptor, pathDirectory.result(name_dataset + "/" + type_pocket) + "global_descriptor_" + type_protomol +"_" + type_pocket, pocket_separation = option_separate )
    return dico_descriptor


def loadDescriptorPocket (dir_pocket="", dico_load = {}, type_descriptor = "", type_protomol = "", type_pocket = "", path_desc_file = "", debug = 0):
    """
    Load descriptors -> load in dictionary
    args: -> path directory pocket in result
          -> structure dico_load
          -> type descriptor (energy, area, volume, ...)
          -> type protomol use (Fpocket or Surflexe)
    return: NONE -> dynamic change directly dico_load
    """
    if path_desc_file == "" : 
        path_file_descriptor = dir_pocket + "descriptor_" + type_descriptor + "_" +type_protomol + "_" + type_pocket
    else : 
        path_file_descriptor = path_desc_file
        
        
    if debug : print os.path.getsize(path_file_descriptor)
    if not os.path.exists(path_file_descriptor) :
        if debug : print "ERROR -> " + path_file_descriptor
        return
    else : 
        filin = open (path_file_descriptor, 'r')
        list_lines = filin.readlines ()
        filin.close ()    
    
    if not type_descriptor in dico_load.keys () : 
        dico_load[type_descriptor] = {}
        dico_load = dico_load[type_descriptor]
    else : 
        dico_load = dico_load[type_descriptor]
     
    for line in list_lines : 
        if re.search ("^pocket",line) : 
            in_dico = line.split ("\t") [0]
#           if debug : print in_dico, dir_pocket
            value = line.split ("\t") [1].strip()
            in_dico = in_dico[7:]
            if in_dico == "pock_vol" : 
                in_dico = "Real_volume"
            if in_dico == "%_ATOM_CONVEXE" : 
                in_dico = "X._ATOM_CONVEXE"
            if in_dico == "CONVEX-SHAPE_COEFFICIENT" : 
                in_dico = "CONVEX.SHAPE_COEFFICIENT"
#            if in_dico == "lig_vol" or in_dico == "pock_vol": #######pass for volume
#                continue
#           if debug : print in_dico
            if not in_dico in dico_load : 
                dico_load[in_dico] = []
            dico_load[in_dico].append (value)

                
def loadDescriptorFPocket (dir_pocket, dico_load, type_protomol, type_pocket, debug = 0):
    """
    A MODIFIER POUR LES CAS DE CHARGEMENT DANS LE FICHIER DE POCHE
    Load descriptors -> load in dictionary
    args: -> path directory pocket in result
          -> structure dico_load
          -> type descriptor (energy, area, volume, ...)
          -> type protomol use (Fpocket or Surflexe)
    return: NONE -> dynamic change directly dico_load
    """
    list_file = os.listdir(dir_pocket)
    for file_descriptor in list_file :
        if re.search ("atm.pdb", file_descriptor) : 
            path_file_fpocket = dir_pocket + file_descriptor
        if re.search ("descriptor_fpocket", file_descriptor) : 
            path_file_descriptor = dir_pocket + file_descriptor

    if "path_file_descriptor" in locals ():
        loadDescriptorPocket (dir_pocket, dico_load, "fpocket", type_protomol, type_pocket)
        
    else : 
        if not "fpocket" in dico_load.keys () : 
            dico_load["fpocket"] = {}
            dico_load = dico_load["fpocket"]
        else : 
            dico_load = dico_load["fpocket"]
        
        
        list_file = os.listdir(dir_pocket)
        filin = open (path_file_fpocket, "r")
        list_lines = filin.readlines ()
        filin.close ()
        list_descriptors = StructureStock.listDescriptorFpocket()
        for line in list_lines : 
            if re.search ("HEADER", line) : 
                for descriptor in list_descriptors : 
                    if re.search (descriptor, line) : 
                        value = line.strip().split (":")[-1].replace(" ", "")
                        if debug : print value, descriptor
                        descriptor_in = descriptor.replace (" ", "_")
                        if not descriptor_in in dico_load.keys () : 
                            dico_load[descriptor_in] = []
                        dico_load[descriptor_in].append (value)



##########def loadDescriptorRadi (dir_pocket, dico_load, type_protomol, type_pocket) : 
##########    """
##########    A reprendre pour pour ajouter PSI/QMD/PCI/INR
##########    Load descriptor from radi3
##########    args: -> PDB ID
##########          -> dictionary descriptor
##########          -> name dataset
##########          -> type of descriptor
##########    return: NONE append radi descriptor in dictionary descriptors
##########    """
##########    
##########    loadDescriptorPocket (dir_pocket, dico_load, "radi", type_protomol, type_pocket)
###########    
###########    path_file_radi = pathDirectory.descriptor() + name_dataset + "/descriptorRadii.txt"
###########    path_file_QMD = pathDirectory.descriptor() + name_dataset + "/radi_QMD"
###########    path_file_PSI = pathDirectory.descriptor() + name_dataset + "/radi_PSI"
###########    path_file_PCI = pathDirectory.descriptor() + name_dataset + "/radi_PCI"
###########    path_file_INR = pathDirectory.descriptor() + name_dataset + "/radi_INR"
###########
###########    # may be modulable path
###########    if not os.path.exists(path_file_radi) : 
###########        print "ERROR ->", path_file_radi, "MISS"
###########        print "Continue"
###########        return
###########    # check key exist
###########    if not type_descriptor in dico_load.keys () : 
###########        dico_load[type_descriptor] = {}
###########    
###########    dico_load = dico_load[type_descriptor]
###########    pocket_result = parseDescriptorRadii3.parseSpecificPocket(PDB_ID, path_file_radi, path_file_QMD ,path_file_PSI, path_file_PCI, path_file_INR) # dictionary Radi
###########    
###########    for descriptor in pocket_result.keys () :
###########        if len (pocket_result[descriptor]) == 1 : # PB case with more data
###########            if not descriptor in dico_load.keys () : 
###########                dico_load[descriptor] = []
###########            dico_load[descriptor] = dico_load[descriptor] + pocket_result[descriptor]
###########            

def loadDescriptorSpecificPDBID (dico_descriptor, list_PDB_ID):
    """
    Write files only list ID PDB
    args: -> dictionary with descriptors
          -> list PDB ID
    return: -> dictionary with only PDB
    """
    
    dico_out = {}
    list_type_pocket = dico_descriptor.keys ()
    list_type_pocket.remove ("data")
    # ionitialization out put
    for type_pocket in list_type_pocket : 
        dico_out[type_pocket] = {}
        list_type_descriptor = dico_descriptor[type_pocket].keys ()
        list_type_descriptor.remove("data")
        for type_descriptor in list_type_descriptor : 
            print type_descriptor
            if type_descriptor == "PDB" : 
                dico_out[type_pocket][type_descriptor] = []
            else :
                dico_out[type_pocket][type_descriptor] = {}
                for descriptor in dico_descriptor[type_pocket][type_descriptor].keys () : 
                    dico_out[type_pocket][type_descriptor][descriptor] = []
    
    # append structure in out put
    for type_pocket in list_type_pocket : 
        for i in xrange (len(dico_descriptor[type_pocket]["PDB"])) : 
            if dico_descriptor[type_pocket]["PDB"][i] in list_PDB_ID : 
                dico_out [type_pocket]["PDB"].append(dico_descriptor[type_pocket]["PDB"][i])
                list_type_descriptors = dico_descriptor[type_pocket].keys () 
                list_type_descriptors.remove ("PDB")
                list_type_descriptors.remove("data")
                for type_descriptor in list_type_descriptors : 
                    for descriptor in dico_descriptor[type_pocket][type_descriptor].keys () : 
                        dico_out [type_pocket][type_descriptor][descriptor].append (dico_descriptor[type_pocket][type_descriptor][descriptor][i])
    return dico_out


def FpocketDruggability (dico_dataset, pocket_retrive_type, name_dataset, only_one_pocket = 1) : 
    """
    
    """
    if name_dataset == "ApoForm128" or name_dataset == "ApoForm138": 
        name_dataset = "ApoForm"
    if only_one_pocket == 0 : 
        print "ERROR -> only one pocket by PDB"
        
    dir_pocket = pathDirectory.descriptor(name_dataset + "/" + pocket_retrive_type)
    
    list_PDB = dico_dataset.keys ()
    for PDB in list_PDB : 
        dir_PDB = dir_pocket + PDB + "/"
        
        list_dir = os.listdir(dir_PDB)
        
        for file_dir in list_dir : 
            if re.search ("^pocket", file_dir) :
                path_dir_pocket = dir_PDB + file_dir 
                path_file_pocket = path_dir_pocket + "/pocket_Fpocket_atom.pdb"
                dico_dataset[PDB]["Drug Score"] = retrieveDrugScore (path_file_pocket)
    
    return dico_dataset
                
                

def retrieveDrugScore (path_file_Fpocket):
    
    filin = open (path_file_Fpocket, "r")
    list_line = filin.readlines ()
    
    for line in list_line : 
        if re.search("HEADER 1", line) : 
            drug_score = line.strip().split (":")[1]
            
    return drug_score




def DogSiteGlobal (p_dir_dataset, l_pdb) : 
    
    dico_out = {}
    for pdb in l_pdb : 
        dico_out[pdb] = {}
#         print pdb
        l_file_pocket = os.listdir (p_dir_dataset + pdb + "/DOGSITE/")
        for p_file in l_file_pocket : 
            if re.search("PocXls_descriptor", p_file) : 
                loadFileDescriptorDogSite (p_dir_dataset + pdb + "/DOGSITE/" + p_file, dico_out[pdb])
                
    return dico_out
                
   
def radiglobal (p_dir_dataset, l_pdb, estimator) : 
    
    dico_out = {}
    for pdb in l_pdb : 
        dico_out[pdb] = {}
        if estimator == "Fpocket" : 
            p_dir_pocket = p_dir_dataset + pdb + "/protein_out/pockets/"
            l_file_pocket = os.listdir (p_dir_pocket)
        else : 
            p_dir_pocket = p_dir_dataset + pdb + "/" + estimator + "/"
            l_file_pocket = os.listdir (p_dir_dataset + pdb + "/" + estimator + "/")
        
        
        for p_file in l_file_pocket : 
            if re.search("_radi_", p_file) : 
                dico_out[pdb][p_dir_pocket + p_file.split ("descriptor_radi_")[0] + ".pdb"] = {}
                loadFileDescriptorRadi (p_dir_pocket + p_file, dico_out[pdb][p_dir_pocket + p_file.split ("descriptor_radi_")[0] + ".pdb"])
                
    return dico_out             
        
        
        
def loadFileDescriptorDogSite (p_filin, dico) : 
    
    filin = open (p_filin, "r")
    l_lines = filin.readlines ()
    filin.close ()
    
    l_desc = l_lines[0].strip().split ("\t")
    
    
        
    for pocket in l_lines[1:] : 
        l_value_desc = pocket.strip().replace (" ", "").split ("\t")
        p_file_pocket = os.path.dirname(p_filin) + "/atms_" + pocket.split ("\t")[0] + ".pdb"
        dico[p_file_pocket] = {}
        
        
        nb_des = len (l_desc)
        i = 0
        while i < nb_des : 
            dico[p_file_pocket][l_desc[i]] = l_value_desc[i]
            i = i + 1 
            
    
            
def loadFileDescriptorRadi (p_filin, dico) : 
    
    filin = open (p_filin, "r")
    l_lines = filin.readlines ()
    filin.close ()
    
    for l in l_lines : 
        desc = l.strip().split ('\t')[0][7:]
        dico[desc] = l.strip().split ('\t')[1]
    
            
         
def FpocketGlobal (p_dir_dataset, l_pdb) : 
    
    dico_out = {}
    for pdb in l_pdb : 
        dico_out[pdb] = {}
        p_dir_pocket = p_dir_dataset + pdb + "/protein_out/pockets/"
        l_file_pocket = os.listdir (p_dir_pocket)
        for p_file in l_file_pocket : 
            if re.search("atm.pdb", p_file) : 
                dico_out[pdb][p_dir_pocket + p_file] = {}
                parseFpocket.resultPocket (p_dir_pocket + p_file, dico_out[pdb][p_dir_pocket + p_file])
    return dico_out           
        
            

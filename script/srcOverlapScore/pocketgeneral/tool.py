import pathDirectory
from os import listdir, system, path
from re import search, sub
from copy import deepcopy
import parsePDB
import writePDBfile

def loadFileDataset (path_file_dataset):
    """
    Generate dictionary dataset with file
    args: -> file dataset
    return: -> dictionary with first key -> PDB ID / second key (ligands and druggability)
    """
    
    dico_dataset = {}
    filin = open (path_file_dataset, "r")
    list_lines = filin.readlines ()
    filin.close ()
    
    for line_dataset in list_lines : 
        PDB_ID = line_dataset.split ("\t")[0].upper ()
        list_ligand = line_dataset.split ("\t")[2].strip().split(" ") # .strip (dell last " ")
        dico_dataset[PDB_ID] = {}
        dico_dataset[PDB_ID]["ligands"] = list_ligand
        dico_dataset[PDB_ID]["druggability"] = line_dataset.split ("\t")[1]
        dico_dataset[PDB_ID]["data"] = line_dataset.strip().split ("\t")[3]
    #print dico_dataset
    return dico_dataset

def openFileforListLines (path_file):
    
    filin = open (path_file, "r")
    list_lines = filin.readlines ()
    filin.close ()
    return list_lines


def listDescriptor (dico_in):
    """
    Generate list descriptor
    args-> dictionary with descriptor
    return -> list descriptors
    """
    
    dico_load = deepcopy(dico_in)
    
    list_descriptor = []
    list_type_descriptor = dico_load.keys ()
    list_type_descriptor.remove ("PDB")
    try : list_type_descriptor.remove ("data")
    except : pass
    for key_dico in list_type_descriptor : 
        for descriptor in dico_load[key_dico] : 
            list_descriptor.append (descriptor)
    
    return list_descriptor


def checkNumberValues (dico_descriptors) : 
    """
    Control function
    test pocket
    """
    
    for key in dico_descriptors.keys () : 
        print key, len (dico_descriptors[key])
    

def generateFilePDBWithChain (dico_dataSet, path_filout, retrieve_type_pocket, name_dataset):
    """
    Generate file with PDB and chain
    args: -> dictionary dataset
          -> path file out 
    return: NONE file out
    """
    
    filout = open (path_filout, "w")
    for PDB_ID in dico_dataSet.keys () :
        filout.write (PDB_ID + "\t") 
        list_chain = []
        list_dir_pocket = pathDirectory.generateListDirPocket(PDB_ID, retrieve_type_pocket, name_dataset)  
        for dir_pocket in list_dir_pocket : 
            list_files = listdir(dir_pocket)
            for file_pocket in list_files : 
                if search (".res",file_pocket) : 
                    file_res = open (dir_pocket + file_pocket,"r")
                    list_lines_res = file_res.readlines ()
                    for line_res in list_lines_res : 
                        chain = line_res.strip()[-1]
                        try:
                            chain = int (chain)
                            continue
                        except :
                            if not chain in list_chain : 
                                filout.write (chain)
                                list_chain.append (chain)
        filout.write ("\n")
    filout.close ()
                        
    
def concatenePocket (list_PDB, path_filout, name_dataset, retrieve_type_pocket, option_only_on_pocket = 1) :   
    """
    Generate one file with every pocket
    args: -> list PDB files
          -> path file out
    return: NONE write file
    """
    
    filout = open (path_filout, "w")
    for PDB_ID in list_PDB : 
        list_dir_pocket = pathDirectory.generateListDirPocket(PDB_ID, retrieve_type_pocket, name_dataset)
        if option_only_on_pocket == 1 and len(list_dir_pocket) > 1 : 
            continue
        else : 
            for pocket_dir in list_dir_pocket : 
                path_file_pocket = pathDirectory.searchPocketAtomFpocket(pocket_dir)
                file_pocket = open(path_file_pocket)
                list_atom_pocket = file_pocket.readlines()
                filout.write("HEADER " + str(PDB_ID) + "\n")
                for atom_pocket in list_atom_pocket : 
                    if search ("^ATOM", atom_pocket) :
                        filout.write (atom_pocket)
                filout.write ("END\n")
    filout.close ()
                    
def fusionDictionaryPocketType (dico_descriptor, data_separated = 1) : 
    """
    Fusion in only one dictionary pocket druggable and no druggable
    args: -> dictionary 
    return: -> dictionary with dictionary
    """
    list_type_pocket = dico_descriptor.keys ()
    
    try : 
        list_type_pocket.remove ("data")
    except : 
        pass
    
    list_PDB = []
    list_data = []
    dico_out = {}
    
    
    
    for type_pocket in list_type_pocket : 
        list_type_descriptor = dico_descriptor[type_pocket].keys ()
#        print list_type_descriptor
        list_PDB = list_PDB +  dico_descriptor[type_pocket]["PDB"]
        list_type_descriptor.remove ("PDB")
        if data_separated : 
            list_data = list_data + dico_descriptor[type_pocket]["data"]
            list_type_descriptor.remove ("data")
            
        for type_descriptors in list_type_descriptor :
            if not type_descriptors in dico_out.keys () :  
                dico_out[type_descriptors] = {}
            for descriptor in dico_descriptor[type_pocket][type_descriptors].keys () : 
                if not descriptor in dico_out[type_descriptors].keys () : 
                    dico_out[type_descriptors][descriptor] = []
                dico_out[type_descriptors][descriptor] = dico_out[type_descriptors][descriptor] + dico_descriptor[type_pocket][type_descriptors][descriptor]
    dico_out["PDB"] = list_PDB
    if data_separated : dico_out["data"] = list_data
    
    return dico_out
                
    
def checkSameOrderDescriptor (dico_descriptors) :
    """
    Check if same descriptor with type of pocket
    args: dictionary descriptor
    return: -> 1 if same descriptor
            -> 0 if not same descriptor
    """
    list_global = []
    list_type_pocket = dico_descriptors.keys ()
    try : list_type_pocket.remove ("data") # remove type data set
    except: pass
    for type_pocket in list_type_pocket : 
        list_descriptor = []
        list_type_descriptor = dico_descriptors[type_pocket].keys ()
        # case empty type
        
        list_type_descriptor.remove ("PDB")
        try : list_type_descriptor.remove ("data")
        except : pass

        if len (list_type_descriptor) == 0 : 
            return 2
        for type_descriptor in list_type_descriptor : 
            for descriptor in dico_descriptors[type_pocket][type_descriptor].keys () : 
                list_descriptor.append (descriptor)
        list_global.append (list_descriptor)
    
    nb_descriptor = len (list_global)
    i=0
    while (i < nb_descriptor - 1) : 
        j = i + 1
        if list_global[i] != list_global[j] : 
            return 0
        i = i + 1
    
    return 1
    
    
def nbValue (dico_descriptor) : 
    """
    Retrieve number of value in dictionary
    arg: -> dictionary with descriptors
    return: -> list of value
            -> 0 if error
    """
    list_nb_value = []
    list_type_descriptor = dico_descriptor.keys ()
    list_type_descriptor.remove ("PDB")
    try : list_type_descriptor.remove ("data")
    except : pass
    for type_descriptor in list_type_descriptor : 
        for descriptor in dico_descriptor[type_descriptor].keys () : 
            nb_value = len (dico_descriptor[type_descriptor][descriptor])
            if not nb_value in list_nb_value : 
                list_nb_value.append (nb_value)
    if len (list_nb_value) == 1 : 
        return list_nb_value[0]
    else : 
        return 0
    

def suppNA (list_in):    
    """
    Dell NA in list values
    args: -> list values
    return: NONE change directly list with values
    """
    
    nb_element = len (list_in)
    i = 0
    while (i < nb_element) : 
        if list_in[i] == "NA" or list_in[i] == "-nan" : 
            del list_in[i]
            nb_element = nb_element - 1
            continue
        i = i + 1
        
       
def appendNAValue (dico_load, debug = 0) : 
    """
    Append NA for have same number of values
    args: dictionary with descriptor
    return: NONE change directly dictionary descriptors
    """
    dico_nb_value = {}
    max_value = 0
    list_type_descriptor =dico_load.keys () 
    list_type_descriptor.remove ("PDB") 
    list_type_descriptor.remove ("data")
    for type_descriptors in list_type_descriptor :
        dico_nb_value[type_descriptors] = {} 
        for descriptor in dico_load[type_descriptors].keys () :
            nb_value = len (dico_load[type_descriptors][descriptor]) 
            dico_nb_value[type_descriptors][descriptor] = nb_value
            if max_value < nb_value :
                if debug : print  type_descriptors,descriptor
                max_value = nb_value
    
    for type_descriptors in dico_nb_value.keys () :
        for descriptor in dico_nb_value[type_descriptors].keys () : 
            if dico_nb_value[type_descriptors][descriptor] != max_value : 
                if debug : print type_descriptors, descriptor,dico_nb_value[type_descriptors][descriptor], max_value
                dico_load[type_descriptors][descriptor].append ("NA")

def fusionFilesHist (path_dir):
    """
    Fusion file distribution
    args: -> path directory
    return: NONE (convect and fusion pdf)
    """
        
    list_file_out = listdir(path_dir)
    list_hist = []
    for file_out in list_file_out : 
        if search ("hist", file_out) : 
            file_pdf = path_dir + "/" +file_out.split (".")[0] + ".pdf"
            list_hist.append(file_pdf)
            cmd_convert = "convert " + path_dir + "/"+ file_out + " " +  file_pdf
            system (cmd_convert)
            system ("rm " +path_dir + "/"+ file_out)
        
    cmd_fusion = "gs -dNOPAUSE -sDEVICE=pdfwrite -sOUTPUTFILE=" + path_dir +"/distribution_descriptors.pdf " + "-dBATCH " + " ".join(list_hist)
    system (cmd_fusion)
    system ("rm " + " ".join(list_hist))

def loadFileWithList (path_file):
    """
    Load file R with index
    args: -> path file
    return: list element
    """
    list_out = []
    file_open = open (path_file, "r")
    list_lines = file_open.readlines ()
    file_open.close ()
    
    for element in list_lines : 
        list_out.append (element.replace ("\"", "").strip ())
    return list_out

def balanceSamePDB (dico_descriptor1,dico_descriptor2) : 
    
    print len (dico_descriptor1["PDB"]), len (dico_descriptor2["PDB"])
    
    list_union1 = list(set(dico_descriptor1["PDB"])-set( dico_descriptor2["PDB"]))
#    print len (list_union1)
    list_union2 = list(set(dico_descriptor2["PDB"])-set( dico_descriptor1["PDB"]))
#    print len (list_union2)
    
    list_intersect = list_union1 + list_union2
#    print list_intersect
    
    i = 0
    while len (dico_descriptor1["PDB"]) != len (dico_descriptor2["PDB"]) : 
        if i < len (dico_descriptor1["PDB"]) : 
            if dico_descriptor1["PDB"][i] in list_intersect :
                suppValue (dico_descriptor1, i)
                continue 
        
        if i < len (dico_descriptor2["PDB"]) : 
            if dico_descriptor2["PDB"][i] in list_intersect :
                suppValue (dico_descriptor2, i)
                continue
        i = i + 1
    
    del dico_descriptor1["PDB"]
    del dico_descriptor2["PDB"]
    
def suppValue (dictionary_descriptor, index, debug = 0):
    
    
    for type_descriptor in dictionary_descriptor.keys () : 
        if type_descriptor == "PDB" : 
            del dictionary_descriptor[type_descriptor][index]
        elif type_descriptor == "data" : 
            del dictionary_descriptor[type_descriptor][index]
            
        else : 
            for descriptor in  dictionary_descriptor[type_descriptor] : 
                if debug : 
                    print descriptor, type_descriptor
                    print len(dictionary_descriptor[type_descriptor][descriptor]), index, type_descriptor, descriptor
                del dictionary_descriptor[type_descriptor][descriptor][index]
                

def selectData (dico_descriptor, IC_descriptor, debug = 0) : 
    
    dico_out = deepcopy(dico_descriptor)
#    print dico_out.keys ()
    
    del dico_out["data"]
    
    for type_pocket in dico_out.keys () : 
        nb_pocket = len (dico_out[type_pocket]["PDB"])
        
        i=0
        while (i < nb_pocket) : 
            in_IC = 0
            for descriptor in IC_descriptor.keys () :
                for type_descriptor in dico_out[type_pocket].keys () :  
                    if descriptor in dico_out[type_pocket][type_descriptor] : 
                        break
                try : 
                    if float(dico_out[type_pocket][type_descriptor][descriptor][i]) > float(IC_descriptor[descriptor]["Borne inf"]) : 
                        if float(dico_out[type_pocket][type_descriptor][descriptor][i]) < float(IC_descriptor[descriptor]["Borne sup"]) : 
                            in_IC = in_IC + 1
                except :
                    pass
            if debug : print in_IC
            if in_IC < int(len (IC_descriptor.keys ())/2) : 
                suppValue (dico_out[type_pocket], i)
                nb_pocket = nb_pocket - 1
            else :
                i = i + 1
        
#        print len (dico_out[type_pocket]["PDB"])
    
    return dico_out
    
            
def generateStructCompositionAtomistic (max_distance, step):
    
    dico_out = {}
    for distance in range (5, max_distance + 1, step) :
        dico_out[distance] = {}
        dico_out[distance]["atom"] = 0
        dico_out[distance]["side_chain"] = 0
        dico_out[distance]["main_chain"] = 0
        dico_out[distance]["sulfur"] = 0 
        dico_out[distance]["carbon"] = 0 
        dico_out[distance]["nitrogen"] = 0
        dico_out[distance]["oxygen"] = 0
        dico_out[distance]["hydrogen"] = 0
        dico_out[distance]["hbond_acceptor"] = 0
        dico_out[distance]["hbond_donor"] = 0
        dico_out[distance]["aromatic"] = 0
        dico_out[distance]["alcool"] = 0
        dico_out[distance]["hydrophobic"] = 0
    
    return dico_out
        
        
    

def transformAA (aa):
    
    aa = aa.upper()
    dico_code = {"S":"SER", "T":"THR", "N":"ASN", "Q":"GLN", "E":"GLU", "D":"ASP", "K":"LYS", "R":"ARG", "H":"HIS", "M":"MET", "C":"CYS", "W":"TRP", "F":"PHE", "Y":"TYR", "A":"ALA", "V":"VAL", "L":"LEU", "I":"ILE", "P":"PRO", "G":"GLY"}
    
    if len (aa) == 1 : 
        return dico_code[aa]
    else :
        for aa_one in dico_code.keys ():
            if dico_code[aa_one] == aa : 
                return  aa_one

def parseSelectDescriptor (path_file_descriptor) : 
    """
    Parse descriptor selected
    arg: -> path file descriptor
    return: -> list descriptor
    """
    
    list_descriptor = []
    filin = open (path_file_descriptor, "r")
    list_line = filin.readlines ()
    filin.close ()
    
    i = 0
    while not search ("accuracy", list_line[i]) :
        descriptor = sub('[ ]{2,}', ' ', list_line[i].strip())
        descriptor = descriptor.replace ("\"", "").split ("]")[-1]
        descriptor = descriptor.split (" ")
        try : 
            descriptor.remove ("")
        except : 
            pass
        
        list_descriptor = list_descriptor + descriptor
        i = i + 1
    
    #acc = [list_line[i+1].split (" ")[1],list_line[i+1].split (" ")[2].strip(),list_line[i+2].split (" ")[-1].strip() ]
    return list_descriptor
    
    
def changeDatasetFileApo (path_dir_descriptor, path_dir_dataset, debug = 1) : 
    
    list_file_descriptor = listdir(path_dir_descriptor)
    
    for file_desc in list_file_descriptor : 
        if search("apo", file_desc) : 
            cmd = "cp " + path_dir_descriptor + file_desc + " " + path_dir_dataset + "protein.pdb"
            
            if debug : print cmd
            system (cmd)
            
    
#def changeDicodatasetHolo (dico_dataset) : 
#    """
#    Change dictionnary dataset
#    arg: -> dico dataset generate current version
#    return: -> dictionnary with only holo forms
#    """
#    
#    dico_out = {}
#    for PDB in dico_dataset.keys () : 
#        if "PDB Holo" in dico_dataset[PDB].keys () : 
#            dico_out [dico_dataset[PDB]["PDB Holo"][0]] = {}
#            dico_out [dico_dataset[PDB]["PDB Holo"][0]] ["druggability"] = dico_dataset[PDB]["druggability"]
#            dico_out [dico_dataset[PDB]["PDB Holo"][0]] ["ligands"] = dico_dataset[PDB]["ligands"]
#            dico_out [dico_dataset[PDB]["PDB Holo"][0]] ["Type structure"] = "holo structure"
#    
#    return dico_out


def selectOnlyTypeStructure (dico_dataset, form_protein) : 
    """
    Change dico dataset, remove holo form
    arg: -> dico dataset generate current version
    return: -> dictionnary with only apo forms
    """
#     print "++++++++++++++"
#     print dico_dataset
#     print "++++++++++++++"
    
    dico_dataset_out = deepcopy(dico_dataset)
    for PDB in dico_dataset_out.keys () : 
#        print dico_dataset_out[PDB]
        if dico_dataset_out[PDB]["Type structure"] != form_protein :
            del dico_dataset_out[PDB]
    
    return dico_dataset_out
    


            
def checkENDFinalLinePDBfile (path_filin_PDB):
    """
    Check if character END in end file PDB
    args: -> path file PDB
    return : -> NONE change directly file PDB
    """
    filin = open (path_filin_PDB, "r")
    list_lines = filin.readlines ()
    filin.close()
    if not search ("^END", list_lines[-1]) :
        if (list_lines[-1][-1] == "\n") :
            filout = open (path_filin_PDB,"a")
            filout.write ("END\n")
            filout.close ()
        else:
            filout = open (path_filin_PDB,"a")
            filout.write ("\nEND\n")
            filout.close ()
    elif not search ("^END\n", list_lines[-1]) :
        filout.write ("\n")
        filout.close ()
    else :
        pass
        
         
def checkHEADERinitialLinePDBfile (path_filin_PDB) : 
    """
    Check if character HEADER in beginning file PDB
    args: -> path file PDB
    return : -> NONE change directly file PDB
    """
    
    filin = open (path_filin_PDB, "r")
    element = filin.read ()
    filin.close()
    if not search ("^HEADER", element) :
        filout = open (path_filin_PDB,"w")
        filout.write ("HEADER \n")
        filout.write (element)
        filout.close ()   
  
  
def formatLinesDescriptor (string_line): # a voir si pas deja la fonction
    
    string_line = sub("[ ]{2,}", " ", string_line)
    string_line = string_line[3:].replace ("\"", "")
    list_descriptor = string_line.split ()
    
    return list_descriptor



def formatLinesAcc (string_line): # a voir si pas deja la fonction
    
    string_line = sub("[ ]{2,}", " ", string_line)
    string_line = string_line.strip()[4:].replace ("\"", "")
    list_acc = string_line.split ("---")
    
    return list_acc


  
def retrieveHoloForm (dico_dataset) : 
    
    list_holo = []
    for PDB_apo in dico_dataset.keys () :
        if not "PDB Holo" in dico_dataset[PDB_apo].keys ():
            continue
        else :
            list_PDB_holo = dico_dataset[PDB_apo]["PDB Holo"]
            for PDB_holo in list_PDB_holo : 
                if not PDB_holo in list_holo : 
                    list_holo.append (PDB_holo)
    
    return list_holo 

def retrieveApoForm (dico_dataset, holo_form):
    
    list_out = []
    for apo_PDB in dico_dataset.keys () :
#        print apo_PDB
        if not "PDB Holo" in dico_dataset[apo_PDB].keys () : 
            pass
        else : 
            if holo_form in dico_dataset[apo_PDB]["PDB Holo"] : 
                if not apo_PDB in list_out : 
                    list_out.append (apo_PDB)
        
    return list_out
            
    
    
def diviseDescriptorHolo (dico_descriptors, dico_dataset, holo_PDB, debug = 0) : 
    

    dico_out = deepcopy (dico_descriptors)
    list_PDB_apo = retrieveApoForm (dico_dataset, holo_PDB)
    if debug : print list_PDB_apo, "list APO"
    
    for type_data in dico_out.keys () : 
        #print type_data
        if not type(dico_out[type_data]) is list :
            i_PDB = 0
            nb_PDB = len (dico_out[type_data]["PDB"])
            if debug : print nb_PDB, "number PDB"
            while i_PDB < nb_PDB and dico_out[type_data]["PDB"] != []: 
                #print dico_out[type_data]["PDB"]
                PDB_ID = dico_out[type_data]["PDB"][i_PDB]
                if debug :print PDB_ID, "PDB testing"
                if not PDB_ID in list_PDB_apo : 
                    del dico_out[type_data]["PDB"][i_PDB]
                    try : del dico_out[type_data]["data"][i_PDB]
                    except : pass
                    nb_PDB = nb_PDB - 1
                    for type_descriptor in dico_out[type_data].keys () : 
                        if type(dico_out[type_data][type_descriptor])is dict : # case data keys
                            for descriptor in dico_out[type_data][type_descriptor].keys () : 
                                del dico_out[type_data][type_descriptor][descriptor][i_PDB]
                        else :
                            pass
                else : 
                    i_PDB = i_PDB + 1 
    
    
    return dico_out
                    

def removeChain (path_protein_PDB) : 
    
    path_directory = path.dirname(path_protein_PDB)
    path_filout = path_directory + "/pralign.pdb"
    
    list_atom = parsePDB.loadCoordSectionPDB(path_protein_PDB)
    
    for atom in list_atom : 
        atom["chainID"] = ""
    
    writePDBfile.coordinateSection(path_filout, list_atom, "ATOM")
    
    return path_filout
    
    
    
    
    
def printDico (dic_in):
    
    l_k1  = dic_in.keys () 
    for k in l_k1 :
        print "---" + str(k) + "---" 
        if type(dic_in[k]) == dict : 
            l_k2 = dic_in[k].keys ()
            for k2 in l_k2 : 
                print "---" + str(k2) + "---" 
                if type(dic_in[k][k2]) == dict : 
                    l_k3 = dic_in[k][k2].keys ()
                    for k3 in l_k3 : 
                        print  "---" + str(k3) + "---" 
                        print dic_in[k][k2][k3]
                else : 
                    print dic_in[k][k2]
        else : 
            print dic_in[k]
            
            
def percentageR(l_desc) : 
    
    i = 0
    nb_desc = len(l_desc)
    while i< nb_desc : 
        if search("^X._", l_desc[i]) : 
            l_desc[i] = l_desc[i].replace ("X.","%")
        i = i + 1
            
            
    
               
            

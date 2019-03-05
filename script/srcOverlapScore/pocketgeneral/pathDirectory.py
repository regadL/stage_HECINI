# global modules
from os import makedirs, listdir, path
from re import search



#globals()["dir_initial"] = "./../../" # don't run with naccess

globals()["dir_initial"] = "/home/borrel/druggabilityProject/"


#globals()["dir_initial"] = "/home/alexandre/intership/"


def dataSet ( dir_in="" ):
    """
    Create dataSet Directory
    args: Directory in dataSet
    return: path
    """
    
    dir = dir_initial + "dataSet/"
    try : makedirs( dir, mode=0777 )
    except : pass
    
    if dir_in != "" : 
        dir_in_dataSet = dir + dir_in + "/"
        try : makedirs( dir_in_dataSet, mode=0777 )
        except : pass
        return dir_in_dataSet
    return dir


def logDir ():
    """
    Create dataSet Directory
    args: Directory log
    return: path
    """
    
    dir = dir_initial + "log/"
    try : makedirs( dir, mode=0777 )
    except : pass
    return dir


def result ( dir_in="" ):
    """
    Create result directory
    args: Directory in dataSet
    return: path
    """
    
    dir = dir_initial + "result/" 
    try : makedirs( dir, mode=0777 )
    except : pass
    if dir_in != "" : 
        dir_in_dataSet = dir + dir_in + "/"
        try : makedirs( dir_in_dataSet, mode=0777 )
        except : pass
        return dir_in_dataSet
    
    return dir


def descriptor ( dir_in=""):
    """
    Create descriptor directory
    args: Directory in dataSet
    return: path
    """
    
    dir = dir_initial + "descriptor/"
    try : makedirs( dir, mode=0777 )
    except : pass
    if dir_in != "" : 
        dir_in_dataSet = dir + dir_in + "/"
        try : makedirs( dir_in_dataSet, mode=0777 )
        except : pass
        return dir_in_dataSet
    return dir

def naccess (path_dir_dataset, PDB_ID):
    
    path_dir = path_dir_dataset + PDB_ID + "/NACCESS/"
    try : makedirs( path_dir, mode=0777 )
    except : pass
    return path_dir
    


def generatePathFileInfoFocket( element_list_PDB, path_dataset ) :
    
    if len( element_list_PDB ) > 1 :
        path_file_PDB = path_dataset + element_list_PDB[0] + "/" + element_list_PDB[0] + "_" + element_list_PDB [1] + "_out/" + element_list_PDB[0] + "_" + element_list_PDB [1] + "_info.txt"
    else : 
        path_file_PDB = path_dataset + element_list_PDB[0] + "/" + element_list_PDB[0] + "_out/" + element_list_PDB[0] + "_info.txt"
    return path_file_PDB



def generatePathFile( element_list_PDB, path_dataset, sufixe_name="" ):
    
    if len( element_list_PDB ) > 1 :
        path_file_PDB = path_dataset + element_list_PDB[0] + "/" + element_list_PDB[0] + "_" + element_list_PDB [1] + sufixe_name
    else : 
        path_file_PDB = path_dataset + element_list_PDB[0] + "/" + element_list_PDB[0] + sufixe_name
    return path_file_PDB



def generateListDirPocket (PDB_ID, retrieve_type_pocket, name_dataset):
    """
    Generate path directory pocket with PDB ID
    args: -> PDB ID
    return: list directories of PDB
    """
    path_descriptor = descriptor(name_dataset + "/" + retrieve_type_pocket + "/" + PDB_ID)
        
 #   print path_descriptor, "Path in search pocket files"
    list_files = listdir(path_descriptor)
    list_dir = []
    for file_dir in list_files : 
        if search ("pocket", file_dir) : 
            list_dir.append (path_descriptor + file_dir + "/")
    return list_dir
    
def searchPocketFile (list_dir_pocket) : 
    """
    Search pocket file (atm and vert)
    args: -> path directory pocket
    return: -> path file atm (pocket file)
            -> path file pqr (protomol file)
    """
    
    list_file_pocket = listdir(list_dir_pocket)
    for file_pocket in list_file_pocket : 
        if search ("atm.pdb",file_pocket) : 
            path_file_atm = list_dir_pocket + file_pocket
        if search ("vert.pqr",file_pocket) : 
            path_file_pqr = list_dir_pocket + file_pocket
    try :return path_file_atm, path_file_pqr
    except : return path_file_atm, ""



def searchPocketResFile(dir_in) :
     
    list_file_pocket = listdir(dir_in)
    for file_pocket in list_file_pocket : 
        if search ("res.pdb",file_pocket) : 
            path_file_res = dir_in + file_pocket
            return path_file_res
        
    return 0





def searchProtomolSurflexe (path_dir_pocket) : 
    
    list_files = listdir(path_dir_pocket)
    for file_result in list_files : 
        if search ("protomol.mol2",file_result) : # run every pocket
            return path_dir_pocket + file_result
    return

def searchPocketAtomFpocket (path_dir_pocket) : 
    
    list_files = listdir(path_dir_pocket)
    
    for file_result in list_files : 
        if search ("_atm.pdb",file_result) : # run every pocket
            if not search ("surflexe", file_result) : 
                return path_dir_pocket + file_result   
    return 
   
   
def searchPocketACC (p_directory):
    l_p_files = listdir(p_directory)
    for p_file in l_p_files : 
        if search("ACC.asa", p_file): 
            return p_directory + p_file
    return None
 
 
def searchFpocketDirResult (path_dir_dataset) :
    
    list_files = listdir(path_dir_dataset)
    for file_dir in list_files : 
        if search ("_out", file_dir) : 
            return path_dir_dataset + file_dir + "/pockets/"
    
    
def listPocketFile (PDB_ID, path_dir_dataset) : 
    
    list_out = []
    path_dir_PDB = path_dir_dataset + PDB_ID
    list_file = listdir(path_dir_PDB)
    for dir_file in list_file : 
        if search ("_out",dir_file) : 
            path_pocket = path_dir_PDB + "/" + dir_file + "/pockets/"
            break
    if not "path_pocket" in locals() :
        return []
    
    list_files_pocket = listdir(path_pocket)
    for file_pocket in list_files_pocket : 
        if search("_atm.pdb",file_pocket) : 
            list_out.append(path_pocket + file_pocket)
    
    return list_out
            
            
def FpocketTest () : 
    """
    Create dataSet Directory
    return: path
    """
    
    dir = dir_initial + "FpocketAnalysis/"
    try : makedirs( dir, mode=0777 )
    except : pass
    
    return dir       

def generatePath (path_directory):
    
    try : makedirs( path_directory, mode=0777 )
    except : pass
    
    return path_directory
    

def searchPocketProximity (path_dir_pocket) : 
    
    list_files = listdir(path_dir_pocket)
    for file_result in list_files : 
        if search ("pocket.pdb",file_result) : # run every pocket
            return path_dir_pocket + file_result   
    return 


def searchLigandPDB (PDB_ID, path_dir_pocket) :
    
    list_files = listdir(path_dir_pocket)
    for file_name in list_files : 
        if len (file_name) <= 7 :
            if search (".pdb", file_name) : 
                return path_dir_pocket + file_name
    
    
def searchModel (path_result) : 
    
    list_out = []
    list_file = listdir(path_result)
    for file_dir in list_file : 
        if search ("quality_predict",file_dir) : 
            if not search (".png",file_dir) and not search (".Rdata",file_dir) : 
                list_out.append (path_result + file_dir)
                
    return list_out
            
        
def searchFileProba (path_directory) : 
    
    list_out = []
    list_files  = listdir(path_directory)
    for name_file in list_files :
        if search("LOO_quality_predict", name_file) and not search(".png", name_file): 
            list_out.append (path_directory + "/" + name_file)
    return list_out
        
        
def searchRMSDFile (path_directoy): 
    
    
    list_file = listdir(path_directoy)
    
    for name_file in list_file : 
        if search ("RMSD", name_file) : 
            return path_directoy + name_file
        
        
        
    
def searchPocketAtomASA(PDB_ID, name_dataset, retrieve_type_pocket) : 
    
    list_dir_pocket = generateListDirPocket(PDB_ID, retrieve_type_pocket, name_dataset)
    if list_dir_pocket == [] : return
    list_files = listdir(list_dir_pocket[0])
    
    
    for file_name in list_files : 
        if search ("ACC.asa", file_name) : 
            return list_dir_pocket[0] + file_name
           
           
def generatePathPocketATOM (PDB_ID, name_dataset, pocket_retrieve_type) : 
    
    path_folder = descriptor (name_dataset) + pocket_retrieve_type + "/" + PDB_ID + "/"
    if not path.isdir (path_folder) :
        path_folder = descriptor (name_dataset) + PDB_ID + "/align_out/"
        if not path.isdir (path_folder) :
            return "NONE"
    for element in listdir (path_folder) : 
        if search ("^pocket", element) : 
            path_pocket = path_folder + element + "/"
            return searchPocketFile(path_pocket)[0]
    return "NONE"



def generatePathMatrix (s_PDB, name_dataset) : 
    
    path_rep = descriptor(name_dataset)
    return path_rep + s_PDB + "/matrix.out"
            
        
def DOGSITE ():
    
    return  dir_initial + "DOGSITE/"
        


def searchDogSiteDirResult (path_dir_with_DogSite_result) : 
    
    
    return listdir(path_dir_with_DogSite_result)
    
        
        
        
        
        
        
        
        
    
    
    
    
    
    



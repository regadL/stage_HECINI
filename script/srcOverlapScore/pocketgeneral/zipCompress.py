import zipfile
from re import search
from os import path, mkdir ,listdir



def createDirectory (path_zip):
    """
    Create directory out with PDB ID
    arg : -> Path zip
    return : -> path directory
    """
    
    path_dir_origin = path.dirname(path_zip)
    s_PDB = path_zip[-8:-4]
    
    path_directory = path_dir_origin + "/" + s_PDB.upper () + "/"
    
    try : mkdir(path_directory)
    except : pass
    
    return path_directory
    

def writeUncompress (content, path_out):
    
    filout = open (path_out, "w")
    filout.write (content)
    filout.close ()
    
    return path_out
    

def unCompressZip (path_zip, path_dir_out):
    
    print path_zip
    
    zfile = zipfile.ZipFile(path_zip, "r")
    
    list_file_zip = zfile.namelist()
    print list_file_zip
    
    for file_in_zip in list_file_zip : 
        if search(".zip", file_in_zip) : 
            file_content = zfile.read (file_in_zip)
            path_zipin = writeUncompress(file_content, path_dir_out + file_in_zip.split ("/")[-1])
            unCompressZip (path_zipin, path_dir_out)
        else : 
            file_content = zfile.read (file_in_zip)
            writeUncompress(file_content, path_dir_out + file_in_zip.split ("/")[-1])
    zfile.close()
    



def unCompressFolder (path_folder):
    
    l_file = listdir(path_folder)
    
    for fil in l_file : 
        if search (".zip", fil) : 
            path_zip = path_folder + fil
            path_dir_out = createDirectory(path_folder + fil)
            unCompressZip(path_zip, path_dir_out)

    else :
        pass



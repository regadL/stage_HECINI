from os import listdir, path
import re


def parseInfoFile (path_info_file):
    
    
    regex_int = re.compile(' ([0-9]+) ')
    regex_float = re.compile('[0-9]+\.[0-9]+')
    file_info = open (path_info_file, "r")
    info_read = file_info.read ()
    file_info.close ()
    list_pocket = info_read.split ("Pocket")[1:]
    
    dico_out = {}
    for pocket in list_pocket :
        list_line = pocket.split ("\n")
        num_pocket = regex_int.findall (list_line[0])[0]
        vol_pocket = regex_float.findall (list_line[7])[0]
        dico_out[num_pocket] = {}
        dico_out[num_pocket]["volume"] = vol_pocket
    
    return dico_out


def volumeNbPocket (list_vol, list_nb_pocket, dir_result):
    
    # check exist file
    if not path.exists(dir_result) : 
        return
    list_files = listdir(dir_result)
    print list_files
    
    for file_fpocket in list_files : 
        if re.search ("_info.txt", file_fpocket) : 
            path_file_info = dir_result + file_fpocket
            dico_result = parseInfoFile(path_file_info)
            for pocket in dico_result.keys () :
                list_vol.append (float (dico_result[pocket]["volume"]))
                list_nb_pocket.append (int (len (dico_result)))

def resultDpocket (path_filin) : 
    """
    Parsing file result dpocket
    args: -> path result
    return: -> dictionary parsed
    """
    
    filin = open (path_filin, "r")
    list_line = filin.readlines ()
    
    dico_parsed = {}
    list_descriptor = re.sub('[ ]{2,}', ' ', list_line[0]).split (" ")
    list_values = re.sub('[ ]{2,}', ' ', list_line[1]).split (" ")
    i = 0
    nb_descriptor = len (list_descriptor)
    while i<nb_descriptor :
        dico_parsed[list_descriptor[i]] = list_values[i]
        i = i +1 
    
    return dico_parsed
        
 
def resultPocket (p_file_pocket, dico_struct): 

        
    filin = open (p_file_pocket, "r")
    l_lines = filin.readlines ()
    filin.close ()
    
    
    for l in l_lines[5:] : 
        if re.search ("HEADER", l) : 
        
        
        
            des = re.sub('[ ]{2,}', ' ', l.split ("- ")[-1].split (":")[0])[:-1]
            val = l.strip ().split (" : ")[-1].replace (" ", "")
            
            dico_struct[des] = val
                
        
        
        
        
    
    
    
    
    
    


   

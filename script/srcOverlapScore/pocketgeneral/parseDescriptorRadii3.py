"""
BORREL Alexandre
05-2012
"""

import re

import pathDirectory


def parseResult (path_filin, debug = 0):
    """
    Parsing complete file -----> not USE
    arg: -> path file in
    return: -> dictionary
    """
    
    list_out = []
    filin = open (path_filin, "r")
    file_read = filin.read ()
    filin.close ()
    
    list_pocket_result = file_read.split (" MOLECULE")
    print len (list_pocket_result)
    for pocket_result in list_pocket_result : 
        pocket_parsed = resultPocketRadi (pocket_result)
        list_out.append (pocket_parsed)
    return list_out

def resultPocketRadi (element_read):
    """
    Retrieve result by pocket
    args: -> element read pocket
    return: dictionary with result
    """
    
    dico_out = {}
    list_lines = element_read.split ("\n")
    
    for line_result in list_lines : 
        if re.search ("MOLECULE", line_result) : 
            nb_atom =  retrieveInt (line_result)[1]
            break
    
    for line_result in list_lines : 
        if re.search("^ PERCENTAGES OF INERTIA", line_result) :
            value_list = retrieveFloat (line_result)
            dico_out ["INERTIA 1"] = [value_list[0]]
            dico_out ["INERTIA 2"] = [value_list[1]]
            dico_out ["INERTIA 3"] = [value_list[2]]
        elif re.search ("^ 3D CONVEX HULL", line_result) : 
            value_line = retrieveInt (line_result)
            dico_out["FACE"] = [value_line[0]]
            #dico_out["EDGES"] = [value_line[1]]
            dico_out["% ATOM CONVEXE"] = [float(value_line[2])/int(nb_atom)]      
        elif re.search ("^ 3D HULL : SURFACE", line_result) : 
            value_line = retrieveFloat (line_result)
            dico_out["SURFACE HULL"] = [value_line[0]]
            dico_out["VOLUME HULL"] = [value_line[1]]
        elif re.search ("^ DIAMETER", line_result) :
            dico_out["DIAMETER HULL"] = retrieveFloat (line_result)
        elif re.search("RADIUS   =", line_result) : 
            dico_out["RADIUS HULL"] = retrieveFloat (line_result)
        elif re.search("^ SMALLEST SIZE :", line_result) : 
            dico_out["SMALLEST SIZE"] = retrieveFloat (line_result)
        elif re.search ("^ HEIGHT =", line_result) : 
            dico_out["RADIUS CYLINDER"] = [retrieveFloat (line_result)[1]]
        elif re.search ("^ GEOMETRICAL CONVEX-SHAPE COEFFICIENT", line_result) : 
            dico_out["CONVEX-SHAPE COEFFICIENT"] = retrieveFloat (line_result)
        elif re.search("^ Radius     :", line_result) : 
            dico_out["RADIUS MIN CYLINDER"] = retrieveFloat (line_result)
        elif re.search("^ Half-height:", line_result) : 
            dico_out["HEIGHT MIN CYLINDER"] = retrieveFloat (line_result)        
            
            
    return dico_out
             
            
            
def retrieveInt (line_result):
    
    list_out = []
    regex_int = re.compile(' ([0-9]+) ')
    result_int = regex_int.findall (line_result)
    if len (result_int) == 1 : 
        return [str (result_int[0])]
    else : 
        for result in result_int : 
            list_out.append (result)
    return list_out


def retrieveFloat (line_result, debug = 0):

    list_out = []
    regex_float = re.compile('[0-9]+\.[0-9]+')
    result_float = regex_float.findall (line_result)
    if debug : print result_float
    
    if len (result_float) == 1 : 
        return [str (result_float[0])]
    else : 
        for result in result_float : 
            list_out.append (result)
    return list_out


def parseSpecificPocket (PDB_ID, path_file_radi, path_file_QMD, path_file_PSI, path_file_PCI, path_file_INR, debug = 0):
    """
    Parse files Radi / QMD / PSI / PCI / INR
    args: -> PDB ID
          -> path QMD file
          -> path PCI file
          -> path INR file
          -> path PSI file 
    return: -> dictionary with descriptors radi
    """
    
    filin = open (path_file_radi, "r")
    file_read = filin.read ()
    filin.close ()    
    
    list_pocket_result = file_read.split (" MOLECULE")
    
    for pocket_result in list_pocket_result : 
        list_lines_pocket = pocket_result.split ("\n")
        for line_pocket in list_lines_pocket : 
            if re.search (PDB_ID, line_pocket) : 
                dico_out = resultPocketRadi (pocket_result)
                pocket_ID = retriveIndexPocket (list_lines_pocket)
                if debug : print pocket_ID
                dico_out["QMD"] = [retrieveValue (path_file_QMD, pocket_ID)]
                dico_out["PSI"] = [retrieveValue (path_file_PSI, pocket_ID)]
                dico_out["PCI"] = [retrieveValue (path_file_PCI, pocket_ID)]
                dico_out["INR"] = [retrieveValue (path_file_INR, pocket_ID)]
                return dico_out


def resultPocketPCI(path_file_PCI, dico_descriptor) : 
    
    filin = open (path_file_PCI, "r")
    list_line = filin.readlines ()
    filin.close ()

    for line_file in list_line : 
        if re.search ("POCKET CONVEXITY  INDEX", line_file) : 
            dico_descriptor["PCI"] = retrieveFloat (line_file)
        elif re.search (" POCKET SPHERICITY INDEX", line_file) : 
            dico_descriptor["PSI"] = retrieveFloat (line_file)
    return dico_descriptor



def retriveIndexPocket (list_line_pocket) :
    
    regex_int = re.compile('[0-9]+')
    return int(regex_int.findall (list_line_pocket[0])[0]) - 1
    
    
def retrieveValue (path_file, pocket_ID, debug = 0) : 
    """
    Retrieve for a line the value
    args: -> Path file
          -> pocket ID
    return: value
    """
    filin = open (path_file, "r")
    list_line_value = filin.readlines ()
    filin.close ()
    
    regex_float = re.compile('[0-9]+\.[0-9]+')
    if debug: print regex_float.findall (list_line_value[pocket_ID])
    return regex_float.findall (list_line_value[pocket_ID]) [0]

def fileResultByPocket (path_filin_radi, path_file_PCI) : 
    
    filin = open (path_filin_radi, "r")
    read_file = filin.read ()
    filin.close ()
    
    dico_descriptor = resultPocketRadi(read_file)
#    print dico_descriptor
    return (resultPocketPCI(path_file_PCI, dico_descriptor))

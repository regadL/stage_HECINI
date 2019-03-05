"""
BORREL Alexandre
04-2012
"""

import parsePDB
import writePDBfile
import runOtherProg

from re import search
from os import path, system



def AllPDBSepareChain (PDB_ID, path_dir_dataset):
    """
    Separe PDB by chain
    arg: -> PDB ID
         -> path directory data set
    return: NONE write one PDB by chain in directory Dataset
    """
    
    # separate chain
    PDB_ID = PDB_ID.upper ()
    path_PDB_file = path_dir_dataset + PDB_ID + "/" + PDB_ID + ".pdb"
    separeByChain(path_PDB_file)
        

    
def AllPDBProtonation (PDB_ID,path_dir_dataset): 
    """
    Protonation only complete protein without ligand and heteroatome
    args : -> PDB ID
           -> path directory dataset
    return: empty -> run PDB2PQR
    """
    
    if runOtherProg.pdb2pqr(path_dir_dataset + PDB_ID + "/protein.pdb") == 1 :
        cleanPQR(path_dir_dataset + PDB_ID + "/protein.pqr")


def separeByChain (path_PDB_file):
    
    file_PDB_parsed = parsePDB.loadCoordSectionPDB(path_PDB_file, section="ATOM")
    
    list_chain = []
    file_open_write = {}
    file_open_write["protein"] = open(path_PDB_file[0:-8] + "protein.pdb", "w")
    for atom_PDB in file_PDB_parsed : 
        chain = atom_PDB["chainID"]
        writePDBfile.coordinateStructure(atom_PDB, "ATOM", file_open_write ["protein"] )
        if not chain in list_chain : 
            list_chain.append (chain)
            file_open_write [chain] = open(path_PDB_file[0:-4] + "_" + chain + ".pdb", "w")
            writePDBfile.coordinateStructure(atom_PDB, "ATOM", file_open_write [chain] )
        else : 
            writePDBfile.coordinateStructure(atom_PDB, "ATOM", file_open_write [chain] )
    
    # close files
    for chain in list_chain : 
        file_open_write[chain].close ()


def cleanPQR( path_file_pqr, debug = 0 ):
    """
    Change residues in pqr
    args: -> path file pqr
          -> debug
    return: NONE write change directly pqr
    """
    
    filin = open (path_file_pqr, "r")
    element_file = filin.read()
    filin.close ()
    
    # same path
    filout = open (path_file_pqr, "w")
    
    # change cys
    element_modif = element_file.replace ("1CBDISU", "CB  CYS")
    element_modif = element_modif.replace ("1SGDISU", "SG  CYS")
    
    list_lines = element_modif.split ("\n")
    
    nb_line = len (list_lines)
    i_line = 0
    while (i_line < nb_line) : 
        if search ("SPP", list_lines[i_line]) or  search ("GLUP", list_lines[i_line]):
            if debug : 
                line_parsed = parsePDB.lineCoords(list_lines[i_line])
                print list_lines[i_line], i_line, path_file_pqr
                print line_parsed["resName"]
        if search ("^ATOM", list_lines[i_line]) : 
            line_parsed = parsePDB.lineCoords(list_lines[i_line])
            #print line_parsed["resName"], i_line, path_file_pqr
            if line_parsed["resName"] == "SPP" : 
                list_lines[i_line] = list_lines[i_line].replace ("ASPP", " ASP")
                if line_parsed["name"] != "HD2" : 
                    filout.write (list_lines[i_line] + "\n")
            elif line_parsed["resName"] == "LUP" : 
                list_lines[i_line] = list_lines[i_line].replace ("GLUP", " GLU")
                if line_parsed["name"] != "HE2" :
                    filout.write (list_lines[i_line] + "\n")
            elif line_parsed["resName"] == "GLU" : 
                if line_parsed["name"] != "HE2" : 
                    filout.write (list_lines[i_line] + "\n")
            elif line_parsed["resName"] == "ASP" : 
                if line_parsed["name"] != "HD2" : 
                    filout.write (list_lines[i_line] + "\n")
            else : 
                filout.write (list_lines [i_line] + "\n")
        else : 
            filout.write (list_lines [i_line] + "\n")
        i_line = i_line + 1     
    filout.close ()



def separeChainFasta (PDB_ID, dir_dataset, debug = 0) : 
    """
    Separe fasta by chain to run water for retrieve identity
    arg: -> PDB id
         -> path directory dataset
    return: NONE write file
    """
    path_file_global = dir_dataset + PDB_ID + "/" + PDB_ID + ".fasta"
    
    if not path.exists(path_file_global) : 
        print "ERROR -> File fasta not exist"
    else : 
        filin = open (path_file_global, "r")
        element_file = filin.read ()
        if debug : print element_file
        filin.close ()
        element_file = element_file.split (">")
        if debug :print element_file
        
        for seq_chain in element_file[1:] : 
            chain = seq_chain[5]
            if debug : print chain, "Chain fasta retrieve"
            path_filout = dir_dataset + PDB_ID + "/" + PDB_ID + "_" + str (chain) + ".fasta"
            filout = open (path_filout,"w")
            filout.write (">" + seq_chain)
            filout.close ()
        

    
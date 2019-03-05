"""
BORREL Alexandre
04-2012
"""

import re

def pasteFastaFileGlobal (list_PDB, path_dataset):
    """
    Paste fasta file
    args: -> list PDB
          -> path data set
    return: NONE (write fasta file and txt with name and chain PDB)
    """
    
    path_filout_fasta = path_dataset + "global_sequence.fasta"
    path_filout_txt = path_dataset + "global_sequence.txt"
    filout_fasta = open (path_filout_fasta, "w")
    filout_txt = open (path_filout_txt, "w")
    
    for PDB in list_PDB : 
        filin = open (path_dataset + PDB + "/" + PDB + ".fasta")
        list_line = filin.readlines ()
        filin.close ()
        for line_PDB in list_line : 
            if re.search ("^>", line_PDB) : 
                element = line_PDB.split ("|") [0]
                element = element.replace (":","_")
                filout_fasta.write ("\n" + element + "\n")
                filout_txt.write (element[1:] + "\n")
            else : 
                filout_fasta.write (line_PDB.strip())
    
    filout_fasta.close ()
    filout_txt.close ()
                 



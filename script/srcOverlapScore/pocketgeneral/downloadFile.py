"""
BORREL Alexandre
04-2012
"""
# global module
from urllib import urlretrieve
from os import system

# personal module
     

def importPDB ( PDB_ID , directory, dir_by_PDB = 1, debug=0 ):
    """Retrieve in http://www.pdb.org PDB file
    args: - list_PDB
          - directory out
    return: - error
            - write in log directory errors download
    NB: make log file best
    """
    
    print PDB_ID, "Download"
    adresseSeq = ( "http://www.pdb.org/pdb/files/%s.pdb" % PDB_ID )
    try:
        path_file_pdb = urlretrieve( adresseSeq )
        if debug : print path_file_pdb
        name_PDB = PDB_ID.upper () + ".pdb"
        if dir_by_PDB : 
            directory_out = directory + PDB_ID + "/"
            system ( "mkdir " + directory_out )
        else : 
            directory_out = directory
        cmd = "mv " + path_file_pdb[0] + " " + directory_out + name_PDB
        if debug : print cmd
        system ( cmd )
        print  str( PDB_ID ) + "-> done"
        return 1
    except:
        print str( PDB_ID ) + "-> ERROR DOWNLOAD PDB file"
        return 0
    
    
def importFasta ( PDB_ID , path_dataset, debug=0 ):
    """Retrieve in http://www.pdb.org PDB file
    args: - list_PDB
          - repertory out
    return: - NULL (print error)
    """
    
    print PDB_ID, "download Fasta"
    adresseSeq = ( "http://www.rcsb.org/pdb/files/fasta.txt?structureIdList=%s.pdb" % PDB_ID )
    try:
        path_file_pdb = urlretrieve( adresseSeq )
        if debug : print path_file_pdb
        name_fasta = PDB_ID.upper () + ".fasta"

        cmd = "mv " + path_file_pdb[0] + " " + path_dataset + PDB_ID + "/" + name_fasta
        if debug : print cmd
        system ( cmd )
        return 1
    except:
        print "ERROR download !!"
        return 0

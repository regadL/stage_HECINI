from os import listdir, path
import descriptor



path_dir = "/home/borrel/pocket_file/"


list_file = listdir(path_dir)


for file_pocket in list_file :
    
    if (path.getsize(path_dir + file_pocket) == 0 ) : 
        continue
    file_PDB = open (path_dir + file_pocket, "r")
    list_line_PDB = file_PDB.readlines()
    file_PDB.close ()
    
    filout = open (path_dir + file_pocket, "w")
    if not list_line_PDB[0] == "HEADER\n" : 
        filout.write ("HEADER\n")
    for line_PDB in list_line_PDB : 
        filout.write (line_PDB.strip() + "                                        \n")
    filout.write("END\n")
    filout.close ()
    
    
    try : descriptor.radi(path_dir + file_pocket, path_dir  + file_pocket[0:4], file_pocket[0:4], ".radi")
    except : print "ERROR !!!!! -> " + file_pocket






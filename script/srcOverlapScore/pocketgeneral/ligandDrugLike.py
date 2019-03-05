
import pathDirectory
import runOtherProg




def retrieveSmileDrugLike (dico_dataset, path_filout, debug = 1):
    
    path_smile_code = path_filout + ".smi"
    path_list_ligand = path_filout + ".list"
    filout_smile = open (path_smile_code, "w")
    filout_ligand = open (path_list_ligand, "w")
    
    
    for PDB_ID in dico_dataset.keys () : 
        print PDB_ID
        
        path_dir_pocket = pathDirectory.generateListDirPocket(PDB_ID, "proximity", "krasowski", dir_pocket=1)
        print path_dir_pocket
        path_pdb_ligand = pathDirectory.searchLigandPDB(PDB_ID, path_dir_pocket[0])
        path_file_smi = runOtherProg.babelPDBtoSMI(path_pdb_ligand)
        code_smile = extractSMILECODE (path_file_smi)
        
        if debug : 
            print path_pdb_ligand
            print path_dir_pocket[0]
            print path_file_smi
        filout_ligand.write (PDB_ID + "\t" + dico_dataset[PDB_ID]["druggability"] + "\n")
        filout_smile.write (code_smile + "\n")
    
    filout_ligand.close ()
    filout_smile.close ()
    
    return [path_file_smi, path_list_ligand]


def extractSMILECODE (path_file_smi, debug = 0):
    
    
    filin = open (path_file_smi, "r")
    list_lines = filin.readlines ()
    filin.close ()
    
    smile_code = list_lines[0].split ("\t")[0]
    if debug : print smile_code, "SMILE CODE"
    
    return smile_code
    


def parseLigandDescriptor (path_file_list_ligand, path_file_descriptor_ligand, path_filout) : 
    
    filout = open (path_filout, "w")
    filout.write ("NameLigand\tDruggable\tfailuresLipinski\tLogP\n")
    
    filin_list = open(path_file_list_ligand, "r")
    list_ligand = filin_list.readlines ()
    filin_list.close ()
    
    filin_descriptor = open (path_file_descriptor_ligand, "r")
    list_molecul_descriptor = filin_descriptor.readlines ()
    filin_descriptor.close ()
    
    
    nb_ligand = len (list_ligand)
    
    for i_ligand in xrange (0,nb_ligand) : 
        
        try : 
            a= list_molecul_descriptor[i_ligand+1].split ("\t")[5]
        except : 
            continue
        
        ligand = list_ligand[i_ligand].strip().split ("\t")
        
        if ligand[1] == "d" and int(list_molecul_descriptor[i_ligand+1].split ("\t")[5]) == 0:
            ok = "ok" 
        elif ligand[1] == "n" and int(list_molecul_descriptor[i_ligand + 1].split ("\t")[5]) > 0 : 
            ok = "ok"
        else :
            ok = "nop"
            
        print i_ligand
        print list_molecul_descriptor[i_ligand]
        try : filout.write (str (ligand[0]) + "\t" + str(ligand[1]) + "\t" + str(list_molecul_descriptor[i_ligand+1].split ("\t")[5]) + "\t" + str (list_molecul_descriptor[i_ligand+1].split ("\t")[1]) + "\t" + str (ok) + "\n")
        except : continue
    
    
    
    
    
       
    


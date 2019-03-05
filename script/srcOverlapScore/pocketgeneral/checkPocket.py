# global module
from os import listdir, system, path, mkdir
from re import search
from shutil import copyfile

# personnal modules
import pathDirectory
import parsePDB
import writePDBfile
import runOtherProg
import dataSet
import tool
import superposeStructure





def selectPocketFpocket (dico_dataSet, path_dir_dataset, name_dataset, pocket_type_retrieve, option_aggregation_pocket=1, file_dir_name = 0 ):
    """
    Select pocket with ligands and Fpocket prediction
    args: 
        -> dictionary dataset with PDB_ID
        -> path directory dataset
        -> option aggregation pocket (check close pocket and identic pocket or same pocket on other chain with 100% identity)
    return: NONE (dell directory pocket in result directory)
    """
    
    for PDB_ID in dico_dataSet.keys () : 
#        print PDB_ID, "Retrieve pocket -> Fpocket"
        for ligand in dico_dataSet[PDB_ID]["ligands"] : 
            selectPocketWithLigandFpocket(PDB_ID, ligand, path_dir_dataset + PDB_ID + "/", name_dataset, file_dir_name = file_dir_name )
        if option_aggregation_pocket : 
            # fusion closed pocket
            checkSamePocket (PDB_ID, pocket_type_retrieve, name_dataset) # pocket close
        checkIdenticPocket (PDB_ID, path_dir_dataset, pocket_type_retrieve, name_dataset)




def selectPocketDogSite (dico_dataSet, path_dir_dataset, name_dataset, pocket_type_retrieve, file_dir_name = 0 ):
    """
    Select pocket with ligands and Fpocket prediction
    args: 
        -> dictionary dataset with PDB_ID
        -> path directory dataset
        -> option aggregation pocket (check close pocket and identic pocket or same pocket on other chain with 100% identity)
    return: NONE (dell directory pocket in result directory)
    """
    
    for PDB_ID in dico_dataSet.keys () : 
#        print PDB_ID, "Retrieve pocket -> Fpocket"
        for ligand in dico_dataSet[PDB_ID]["ligands"] : 
            selectPocketWithLigandDogSite(PDB_ID, ligand, path_dir_dataset + PDB_ID + "/DOGSITE/", name_dataset )





def selectPocketFpocketApo (dico_dataset_apo, dico_dataset_holo, dir_dataset, name_dataset) :
     
    dir_descriptor = pathDirectory.descriptor(name_dataset) 
    
    for PDB_ID in dico_dataset_apo.keys () :
        print dico_dataset_apo[PDB_ID].keys ()
        print PDB_ID, dico_dataset_apo[PDB_ID]["PDB holo"]
        
        for ligand in dico_dataset_holo[dico_dataset_apo[PDB_ID]["PDB holo"][0]]["ligands"] : 
            print ligand, "LIGAND"
            selectPocketWithLigandApo(PDB_ID, dico_dataset_apo[PDB_ID]["PDB holo"][0], ligand, name_dataset)
        
        
        
        
        
        
        
        
#         if dico_dataSet[PDB_ID]["Type structure"] == "apo structure" : 
#             path_file_ligand = pathDirectory.searchLigandPDB(PDB_ID, dir_descriptor + PDB_ID + "/")
#             path_dir_pocket = pathDirectory.searchFpocketDirResult(dir_descriptor + PDB_ID + "/")
#             
#             list_file_pocket = listdir(path_dir_pocket)
#             for file_pocket in list_file_pocket : 
#                 if search("_vert.pqr", file_pocket) : 
#                     path_file_protomol = path_dir_pocket + file_pocket
#                     print path_file_protomol
#                     print path_file_ligand
#                     pocketWithLigandFile (path_file_ligand, path_file_protomol, dir_descriptor + PDB_ID + "/")
#                 else : 
#                     pass
    

def retrieveLigandFromHoloForm (dictionary_dataset, name_dataset) : 
    """
    Retrieve ligand in holo protein
    args: - dictionary dataset
          - name dataset
    return: NONE write ligand file in directory
    """
    
    dir_descriptor = pathDirectory.descriptor(name_dataset)
    dir_dataset = pathDirectory.dataSet(name_dataset)
    
    for PDB_ID in dictionary_dataset.keys ():
        if "PDB holo" in dictionary_dataset[PDB_ID].keys () :
            print "--------" 
            print PDB_ID, 
            print dictionary_dataset[PDB_ID].keys ()
            print dictionary_dataset[PDB_ID]["Type structure"]
            print dictionary_dataset[PDB_ID]["PDB holo"][0]
            print "-----"
            
            path_PDB_holo = dir_dataset + str(dictionary_dataset[PDB_ID]["PDB holo"][0]) + "/" + str(dictionary_dataset[PDB_ID]["PDB holo"][0]) + ".pdb"
            print path_PDB_holo
            
            PDB_holo_parsed = parsePDB.loadCoordSectionPDB(path_PDB_holo)
            
            list_ligand_parsed = parsePDB.retrieveLigand(PDB_holo_parsed, dictionary_dataset[PDB_ID]["ligands"][0])

            if len (list_ligand_parsed) > 1 : 
                print "SEVERAL LIGAND -> MANUAL CHECK ", dictionary_dataset[PDB_ID]["PDB holo"][0]
                writePDBfile.coordinateSection(dir_descriptor + PDB_ID + "/" +  dictionary_dataset[PDB_ID]["ligands"][0] + ".pdb", list_ligand_parsed[0], "HETATM")
            else :          
                writePDBfile.coordinateSection(dir_descriptor + PDB_ID + "/" +  dictionary_dataset[PDB_ID]["ligands"][0] + ".pdb", list_ligand_parsed[0], "HETATM")
    
    
    



def selectPocketProximity (dico_dataSet, path_dir_dataset, name_dataset, pocket_type_retrieve, dir_pocket, debug = 0):
    """
    Retrieve pocket, pocket is proximity of ligands
    args: -> dictionary dataset
          -> path directory dataset
          -> name dataset
    return: NONE write dataset
    """

    for PDB_ID in dico_dataSet.keys () : 
        path_directory = pathDirectory.descriptor(name_dataset + "/proximity/" + PDB_ID)
#        print PDB_ID, "Retrieve pocket -> Proximity"
        for ligand in dico_dataSet[PDB_ID]["ligands"] :
            # path complexe protonated
            path_file_complexe = path_dir_dataset + PDB_ID + "/protein.pdb"
            path_file_PDB = path_dir_dataset+ PDB_ID + "/" + PDB_ID + ".pdb"
            
            if debug : print path_file_complexe, ": path file pdb"
            if not path.exists(path_file_complexe) : 
                continue
            selectPocketWithLigandProximity(path_file_complexe, path_file_PDB, ligand, path_directory, PDB_ID)
            
        checkIdenticPocket (PDB_ID, path_dir_dataset, pocket_type_retrieve, name_dataset, dir_pocket)
    
    


def selectPocketWithLigandProximity(path_file_complexe, path_file_PDB, ligand, path_directory, PDB_ID, threshold = 4, debug = 0 ) :
    """
    Generate file pocket from 4A the ligand
    """
    # parse PDB file (perhaps change with PDB parser)
    PDB_parsed = parsePDB.loadCoordSectionPDB(path_file_PDB)
    complexe_parsed = parsePDB.loadCoordSectionPDB(path_file_complexe)
    list_atom_ligands = parsePDB.retrieveLigand(PDB_parsed, ligand)
    
    nb_pocket = 1
    for ligand_parsed in list_atom_ligands :
        dir_pocket = path_directory + "pocket" + str (nb_pocket) + "/"
        pathDirectory.generatePath(dir_pocket) # generate path not return
        nb_pocket = nb_pocket + 1
        generateFilePocket (complexe_parsed, ligand_parsed, dir_pocket, threshold, PDB_ID)
        writePDBfile.coordinateSection(dir_pocket + ligand_parsed[0]["resName"] + ".pdb", ligand_parsed, "HETATM", header= PDB_ID, connect_matrix = 0)
  

def generateFilePocket (complexe_parsed, list_atom_ligand, dir_pocket, threshold, PDB_ID):
    """
    Generate file pocket with PDB parsed
    """
    list_serial = []
    filout = open (dir_pocket + "pocket_atm.pdb", "w")
    filout.write ("HEADER " + PDB_ID + "\n")
    for atom_ligand in list_atom_ligand : 
        for atom_complex in complexe_parsed :
            if not atom_complex["serial"] in list_serial :
                if parsePDB.distanceTwoatoms(atom_ligand, atom_complex) <= threshold :
                    writePDBfile.coordinateStructure(atom_complex, "ATOM", filout)
                    list_serial.append (atom_complex["serial"])
    filout.write ("END\n")
    filout.close ()



def checkIdenticPocket (PDB_ID, path_dir_dataset, retrieve_type_pocket, name_dataset, debug = 0):
    """
    Check identical pocket on different chains
    args: - PDB ID
          - path directory dataset
    return: NONE (dell directly directory in result folder)
    """
    
    list_dir_pocket = pathDirectory.generateListDirPocket(PDB_ID, retrieve_type_pocket, name_dataset)
    if list_dir_pocket == [] : 
        return
    nb_pocket = len (list_dir_pocket)
    if nb_pocket == 1 : 
        return # no check if is only one pocket
    
    i_pocket = 0
    print nb_pocket
    while i_pocket <  nb_pocket :
        path_pocket_atm1, path_pocket_vert1 = pathDirectory.searchPocketFile (list_dir_pocket [i_pocket])
        pocket_parsed1 = parsePDB.loadCoordSectionPDB(path_pocket_atm1)
        i_pocket2 = i_pocket + 1
        while i_pocket2 < nb_pocket :
            path_pocket_atm2, path_pocket_vert2 = pathDirectory.searchPocketFile (list_dir_pocket [i_pocket2])
            pocket_parsed2 = parsePDB.loadCoordSectionPDB(path_pocket_atm2)
            if identicPocket (pocket_parsed1, pocket_parsed2, path_dir_dataset, PDB_ID) == 1 :
                delPocket(list_dir_pocket [i_pocket2] )
                del list_dir_pocket [i_pocket2]
                nb_pocket = nb_pocket - 1
            else :
                i_pocket2 = i_pocket2 + 1
        i_pocket = i_pocket + 1
   
        
    if debug : print list_dir_pocket, "debug"



def identicPocket (pocket1_parsed, pocket2_parsed, path_directory_dataset, PDB_ID, debug = 0):
    """
    Check pockets are identic one 2 chains with 100% identity
    args:
        -> pocket parsed 1 (parsed with parsePDB)
        -> pocket parsed 2 (parsed with parsePDB)
        -> path directory dataset
        -> PDB ID
    return :
        -> 0 no identic pocket
        -> 1 identic pocket
    """
    # parse file identity
#    dico_identity = dataSet.selectComplex2calculDescriptor(path_directory_dataset + "redondance_dataset")
#    if dico_identity == {} : 
    dico_identity = dataSet.calculIdenticWater(PDB_ID, path_directory_dataset)
    
    if debug : print dico_identity
    list_chain_pocket1 = []
    list_chain_pocket2 = []
    for atom1 in pocket1_parsed :
        if not atom1["chainID"] in list_chain_pocket1 : 
            list_chain_pocket1.append (atom1["chainID"])
    for atom2 in pocket2_parsed :
        if not atom2["chainID"] in list_chain_pocket2 : 
            list_chain_pocket2.append (atom2["chainID"])
    if debug : 
        print list_chain_pocket1, "list chain 1"
        print list_chain_pocket2, "list chain 2"
    
#    if len(list_chain_pocket1) != len (list_chain_pocket2) : 
#        return 0

    flag_identic = 0
    for chain1 in list_chain_pocket1 : 
        for chain2 in list_chain_pocket2 : 
            key_chain1 = PDB_ID + "_" + chain1
            key_chain2 = PDB_ID + "_" + chain2
            if debug : 
                print key_chain1, key_chain2, "Check chain tested"
                
            if key_chain1 != key_chain2 :
                if key_chain1 in dico_identity and key_chain2 in dico_identity[key_chain1] :
                    if debug : 
                        print dico_identity[key_chain1][key_chain2]
                    if dico_identity[key_chain1][key_chain2] == "100.0%" :
                        flag_identic = flag_identic + 1
                else : 
                    if dico_identity[key_chain2][key_chain1] == "100.0%" :
                        flag_identic = flag_identic + 1
    
    if debug : print flag_identic                  
    if flag_identic >= 1  : #len (list_chain_pocket1) : # found equivalent chain maybe change
        return 1
    else :
        return 0
                
   
def checkSamePocket (PDB_ID, retrieve_type_pocket, name_dataset, debug = 0):
    """check close pocket and fusion pcket closed
    args: -> PDB ID
    return: 0 or 1"""
    
    list_dir_pocket = pathDirectory.generateListDirPocket(PDB_ID, retrieve_type_pocket, name_dataset)
    nb_pocket = len (list_dir_pocket)
    if nb_pocket == 1 : 
        return
    
    i_pocket = 0
    while i_pocket <  nb_pocket :
        path_pocket1_atm, path_pocket1_pqr = pathDirectory.searchPocketFile (list_dir_pocket [i_pocket])
        pocket_parsed1 = parsePDB.loadCoordSectionPDB(path_pocket1_atm)
        i_pocket2 = i_pocket + 1
        while i_pocket2 < nb_pocket : 
            path_pocket2_atm, path_pocket2_pqr = pathDirectory.searchPocketFile (list_dir_pocket [i_pocket2])
            pocket_parsed2 = parsePDB.loadCoordSectionPDB(path_pocket2_atm)

            if samePocket (pocket_parsed1, pocket_parsed2) == 1 :
                fusionPocket (list_dir_pocket [i_pocket], list_dir_pocket [i_pocket2] )
                del list_dir_pocket [i_pocket2]
                nb_pocket = nb_pocket - 1
            
            i_pocket2 = i_pocket2 + 1
        i_pocket = i_pocket + 1


def minimalDistancePocket (pocket_parsed1, pocket_parsed2)  : 
    
    distance = 1000
    for atom1 in pocket_parsed1 : 
        for atom2 in pocket_parsed2 : 
            distance_temp = parsePDB.distanceTwoatoms(atom1, atom2)
            if distance_temp < distance :
                distance = distance_temp
    return float (distance) 
    
            
def samePocket (pocket1, pocket2, distance = 4):
    
    for atom1 in pocket1 : 
        for atom2 in pocket2 : 
            if parsePDB.distanceTwoatoms(atom1, atom2) < distance :
                return 1
    return 0

def fusionPocket (path_dir_pocket1, path_dir_pocket2, debug = 0):
    
    # select file
    path_pocket1_atm, path_pocket1_pqr = pathDirectory.searchPocketFile (path_dir_pocket1)
    path_pocket2_atm, path_pocket2_pqr = pathDirectory.searchPocketFile (path_dir_pocket2)
    
    # fusion file (sum value header, and mean)
    fusionFilePocket (path_pocket1_atm, path_pocket2_atm)
    fusionFilePocket (path_pocket1_pqr, path_pocket2_pqr)
    
    cmd = "rm -r " + path_dir_pocket2
    if debug : print cmd
    system (cmd)
    if debug : print "[VERBOSE] -> FUSION" + path_dir_pocket1 + "+" + path_dir_pocket2
    
    
def delPocket (path_dir_pocket, debug = 1):
    
    cmd = "rm -r " + path_dir_pocket
    if debug : print cmd
    system (cmd)   
    
        
def fusionFilePocket (path_pocket1, path_pocket2, debug = 0) :
    """
    Fusion file pocket (merge header)
    args: -> path file pocket 1
          -> path file pocket 2
    return: NONE (remove file pocket 2 and append element in file pocket 2)
    """
    
    if debug : print "FUSION -> ", path_pocket1, path_pocket2 
    filin1 = open (path_pocket1, "r")
    list_lines1 = filin1.readlines ()
    filin1.close ()
    filin2 = open (path_pocket2, "r")
    list_lines2 = filin2.readlines ()
    filin2.close ()
    filout = open (path_pocket1, "w")
    
    i_line2 = 0
    nb_line2 = len (list_lines2)
    while i_line2 < nb_line2 :
        if search ("HEADER 0", list_lines2[i_line2]) or  search ("HEADER 1", list_lines2[i_line2]) or search ("HEADER 6", list_lines2[i_line2]) or search ("HEADER 7", list_lines2[i_line2]) or search ("HEADER 8", list_lines2[i_line2]) or search ("HEADER 10", list_lines2[i_line2]) or search ("HEADER 11", list_lines2[i_line2]):
            filout.write (list_lines1[i_line2].split (":")[0] + ": NA\n") 
        elif search ("HEADER 2", list_lines2[i_line2]) or  search ("HEADER 9", list_lines2[i_line2]) or search ("HEADER 12", list_lines2[i_line2]) :
            value = float(list_lines2[i_line2].strip().split (":")[-1].replace (" ", "")) + float(list_lines1[i_line2].strip().split (":")[-1].replace (" ", "")) 
            filout.write (list_lines1[i_line2].split (":")[0] + ": " + str (value) + "\n") 
        elif search ("HEADER 3", list_lines2[i_line2]) or search ("HEADER 4", list_lines2[i_line2]) or search ("HEADER 5", list_lines2[i_line2]): 
            value = (float(list_lines2[i_line2].strip().split (":")[-1].replace (" ", "")) + float(list_lines1[i_line2].strip().split (":")[-1].replace (" ", "")) ) / 2
            filout.write (list_lines1[i_line2].split (":")[0] + ": " + str (value) + "\n")         
        elif search ("HEADER", list_lines2[i_line2]) : 
            filout.write (list_lines2[i_line2]) 
        elif search("ATOM",list_lines2[i_line2] ) :
            filout.write (list_lines2[i_line2])
        i_line2 = i_line2 + 1
    for line_filin1 in list_lines1 : 
        if not search("HEADER", line_filin1) : 
            filout.write (line_filin1)
    filout.close ()
            

def selectPocketWithLigandFpocket (PDB_ID, name_ligand, path_dir_with_Fpocket_result, name_dataset, file_dir_name = 0, clean = 1, PDB_holo = 0, debug = 0) :
    """
    Select pocket by PDB ID and name ligand, apo and holo conception also
    args: -> PDB ID
          -> ligand ID
          -> path directory
    return: NONE copy files 
    """
    
    # path directory Fpocket used
    path_dir_Fpocket = pathDirectory.searchFpocketDirResult (path_dir_with_Fpocket_result)
    # path directory result
    if file_dir_name != 0 :
        path_dir_descriptor = pathDirectory.descriptor(name_dataset + "/Fpocket/" + PDB_ID)
    else : 
        path_dir_descriptor = pathDirectory.descriptor(name_dataset + "/" + PDB_ID)
    
    if clean : 
        system ("rm -r " + path_dir_descriptor + "*")
    
    if PDB_holo != 0 : 
        path_file_originalPDB = pathDirectory.dataSet(name_dataset + "/" + PDB_holo) + PDB_holo + ".pdb"
    else : 
        path_file_originalPDB = path_dir_with_Fpocket_result + PDB_ID.upper () + ".pdb"
    
    # load atom ligands
    originalPDB_parsed = parsePDB.loadCoordSectionPDB(path_file_originalPDB)
    ligand_parsed = parsePDB.retrieveLigand(originalPDB_parsed, name_ligand)
    index = 1
    for ligand in ligand_parsed :
        writePDBfile.coordinateSection(path_dir_descriptor + name_ligand + "_" + str (index) + ".pdb", ligand, "HETAM")
    
    # search in .pqr (protomol determined by Fpocket) if distance with ligand < 1
    list_files_pocket = listdir(path_dir_Fpocket)
    for file_in_dir_Fpocket in list_files_pocket :
        if file_in_dir_Fpocket [-8:] == "vert.pqr" :
            path_file_protomol = path_dir_Fpocket + file_in_dir_Fpocket
            if debug :
                print path_file_protomol[0:-8]
            protomol_parsed = parsePDB.loadCoordSectionPDB (path_file_protomol)
            
            if ligandConfusedProtomol(protomol_parsed, ligand_parsed) == 1 : 
                path_dir_out = path_dir_descriptor + file_in_dir_Fpocket.split ("_")[0] + "/"
                if debug : 
                    print path_dir_out
                if not path.exists(path_dir_out) : 
                    mkdir(path_dir_out)
                
                cmd = "cp " + path_file_protomol [0:-8] + "* " + path_dir_out
                if debug : 
                    print cmd, "--> cp file pocket in result directory <--"
                system (cmd)





def selectPocketWithLigandDogSite (PDB_ID, name_ligand, p_DogSite_pocket, name_dataset, clean = 1, debug = 0) :
    """
    Select pocket by PDB ID and name ligand, apo and holo conception also
    args: -> PDB ID
          -> ligand ID
          -> path directory
    return: NONE copy files 
    """
    
    # path directory Fpocket used
    path_dir_descriptor = pathDirectory.descriptor(name_dataset + "/DogSite/" + PDB_ID + "/")
    
    if clean : 
        system ("rm -r " + path_dir_descriptor + "*")
    
    path_file_originalPDB = pathDirectory.dataSet(name_dataset + "/" + PDB_ID) + PDB_ID + ".pdb"
    
    # load atom ligands
    originalPDB_parsed = parsePDB.loadCoordSectionPDB(path_file_originalPDB)
    ligand_parsed = parsePDB.retrieveLigand(originalPDB_parsed, name_ligand)
    
    index = 1
    for ligand in ligand_parsed :
        writePDBfile.coordinateSection(path_dir_descriptor + name_ligand + "_" + str (index) + ".pdb", ligand, "HETAM")
    
        # search in pocket file distance ligand < 3
        l_p_pocket = listdir(p_DogSite_pocket)
        i_pocket = 1
        for p_pocket in l_p_pocket :
            
                pocket_parsed = parsePDB.loadCoordSectionPDB (p_DogSite_pocket + p_pocket)
                
                if ligandConfusedProtomol(pocket_parsed, ligand_parsed, 3) == 1 : 
                    
                    p_dir_desc_pocket = path_dir_descriptor + "pocket" + str(i_pocket) + "/"
                    pathDirectory.generatePath(p_dir_desc_pocket)
                    i_pocket = i_pocket + 1 
                    
                    p_file_out = p_dir_desc_pocket + "pocket_DogSite_atm.pdb"
                    
                    
                    # write again the PDB file pocket because class PDB does not run whit this data
                    
                    writePDBfile.coordinateSection(p_file_out, pocket_parsed, recorder = "ATOM")
                    
                    
                    
                    
#                     cmd = "cp " + p_DogSite_pocket + p_pocket + " " + p_file_out
                   
                    
                    
#                     if debug : 
#                         print cmd, "--> cp file pocket in result directory <--"
#                     system (cmd)




def selectPocketWithLigandApo (PDB_apo, PDB_holo, ligand, name_dataset, debug = 1) :
    """
    Select pocket by PDB ID and name ligand, apo and holo conception also
    args: -> PDB ID
          -> ligand ID
          -> path directory
    return: NONE copy files 
    """
    
    # path directory Fpocket used
    path_dir_Fpocket = pathDirectory.searchFpocketDirResult (pathDirectory.dataSet(name_dataset) + PDB_apo + "/")
    print path_dir_Fpocket, "CHECK 443 !!!!"
    
    # path directory result
    path_dir_descriptor = pathDirectory.descriptor(name_dataset + "/Fpocket/" + PDB_apo)
    print path_dir_descriptor, "CHECK 447 !!!!"
    
    # path rotation matrix
    path_matrix = path_dir_descriptor + "matrix.out"
    d_matrix = superposeStructure.formatMatrix(path_matrix)
    
    # retrieve ligand
    path_holo_PDB = pathDirectory.dataSet(name_dataset + "/" + PDB_holo) + PDB_holo + ".pdb"
    holo_parsed = parsePDB.loadCoordSectionPDB(path_holo_PDB)
    l_ligand_parsed = parsePDB.retrieveLigand(holo_parsed, ligand)
    
    for i in xrange (0,len(l_ligand_parsed)) : 
        superposeStructure.applyMatrixLigand (l_ligand_parsed[i], d_matrix)
    
        writePDBfile.coordinateSection(path_dir_descriptor + "ligand_transloc" + str(i) + ".pdb", l_ligand_parsed[i], "HETATM")
    
    
    # search in .pqr (protomol determined by Fpocket) if distance with ligand < 1
    list_files_pocket = listdir(path_dir_Fpocket)
    for file_in_dir_Fpocket in list_files_pocket :
        if file_in_dir_Fpocket [-8:] == "vert.pqr" :
            path_file_protomol = path_dir_Fpocket + file_in_dir_Fpocket
            
            if debug :
                print path_file_protomol[0:-8]
            protomol_parsed = parsePDB.loadCoordSectionPDB (path_file_protomol)
             
            if ligandConfusedProtomol(protomol_parsed, l_ligand_parsed, distance=3) == 1 : # change thresold for apo -> flexibility
                path_dir_out = path_dir_descriptor + file_in_dir_Fpocket.split ("_")[0] + "/"
                if debug : 
                    print path_dir_out, "RESULT"
                if not path.exists(path_dir_out) : 
                    mkdir(path_dir_out)
                 
                cmd = "cp " + path_file_protomol [0:-8] + "* " + path_dir_out
                if debug : 
                    print cmd, "--> cp file pocket in result directory <--"
                system (cmd)





def pocketWithLigandFile (path_file_ligand, path_file_protomol, path_directory, debug = 1):
    
    protomol_parsed = parsePDB.loadCoordSectionPDB(path_file_protomol)
    ligand_parsed = parsePDB.loadCoordSectionPDB(path_file_ligand)
    
    
    if ligandConfusedProtomol(protomol_parsed, [ligand_parsed]) == 1 :
        print  path_file_protomol.split ("_")[-2].split ("/")[-1]
        path_dir_out = path_directory + path_file_protomol.split ("_")[-2].split ("/")[-1] + "/"
        if debug : 
            print path_dir_out, "path directory out pocket selected"
        if not path.exists(path_dir_out) : 
            mkdir(path_dir_out)
            
            cmd = "cp " + path_file_protomol [0:-8] + "* " + path_dir_out
            if debug : 
                print cmd, "--> cp file pocket in result directory <--"
            system (cmd)
    
    
    
    
    
def ligandConfusedProtomol (protomol_parsed, list_ligand_parsed, distance = 1.0):
    """
    Check if ligand is confused with protomol
    args: -> protomol parsed
          -> ligand parsed
    return: -> 0 if not confuse or 1 is confuse
    NB: change with parser PDB.py to limit modules
    """
    
    for atom_protomol in protomol_parsed :
        for ligand_parsed in list_ligand_parsed : 
            for atom_ligand in ligand_parsed : 
                if parsePDB.distanceTwoatoms(atom_ligand, atom_protomol) < distance  :
                    return 1
    return 0
                        

def retrieveListRediduesPocket(path_file_pocket, debug = 0) : 
    """
    Retrieve list residues by pocket for Surflexe
    args: -> path file pocket (.pqr in result pocket folder)
    return: NONE path file generated (.res)
    
    """
    pocket_parsed = parsePDB.loadCoordSectionPDB(path_file_pocket)
    path_filout = path_file_pocket[0:-8] + ".res"
    
    filout = open(path_filout, "w")
    list_res = []
    for atom_pocket in pocket_parsed :
        write_element = atom_pocket["resName"] + str (atom_pocket["resSeq"]) + atom_pocket["chainID"]
        if not write_element in list_res : 
            filout.write (write_element + "\n")
            list_res.append (write_element)
    filout.close ()
    return path_filout


def distancePocket (list_PDB, path_dir_pocket, path_filout_distance):
    """
    Retrieve by PDB and pockets the minimal distance
    args: -> list PDB
          -> path filout_all distance
    return: -> NONE (write file with distance)
    """
    
    filout_all = open (path_filout_distance, "a")
    filout_min = open (path_filout_distance + "_min", "a")
    
    for PDB_ID in list_PDB : 
        list_file_pocket = pathDirectory.listPocketFile (PDB_ID, path_dir_pocket)
        i_pocket1 = 0
        nb_pocket = len (list_file_pocket)
        while (i_pocket1 < nb_pocket - 1) : 
            i_pocket2 = i_pocket1 + 1
            pocket_parsed1 = parsePDB.loadCoordSectionPDB(list_file_pocket[i_pocket1])
            pocket_parsed2 = parsePDB.loadCoordSectionPDB(list_file_pocket[i_pocket2])
            dist_pocket = minimalDistancePocket (pocket_parsed1, pocket_parsed2)
            if dist_pocket < 21 : 
                filout_min.write (str(dist_pocket) + "\n")
                filout_all.write (str(dist_pocket) + "\n")
            else : 
                filout_all.write (str (dist_pocket) + "\n")
                
            i_pocket1 = i_pocket1 + 1
    filout_all.close()
    filout_min.close ()
    
    runOtherProg.runRscriptHisto(path_filout_distance, "minimal_distance", 1)    
    runOtherProg.runRscriptHisto(path_filout_distance + "_min", "minimal_distance", 1) 



def retrievePocketDogSite ( p_dataset, l_pdb, debug = 1):
    
    p_dogsite = pathDirectory.DOGSITE()
    l_dir_pocket = []
    
    for pdb_ID in l_pdb : 
        print pdb_ID
        p_dir_dataSet_pocket = p_dataset + pdb_ID + "/DOGSITE/"
        
        # create folder for pocket
        pathDirectory.generatePath(p_dir_dataSet_pocket)
        
        # move pocket file in dataset directory
        l_file_dogsite = listdir(p_dogsite + pdb_ID + "/")
        for name_file in l_file_dogsite : 
            p_file = p_dogsite + pdb_ID + "/" + name_file
            if search(".pdb", p_file) : 
                copyfile (p_file, p_dir_dataSet_pocket + name_file )
            elif search("PocXls_descriptor", p_file) : 
                if not search (".txt2",p_file ) : 
                    copyfile (p_file, p_dir_dataSet_pocket + name_file )
                    l_dir_pocket.append (p_dir_dataSet_pocket)
    
    return l_dir_pocket
        
        
    
    
    
    
    



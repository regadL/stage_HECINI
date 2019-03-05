from math import sqrt
import os
from PDB import *

import runOtherProg
import pathDirectory
import parseNACCESS



def globalGetCompositionPocket (dico_dataset, path_dir_dataset, type_pocket, name_dataset, debug = 0):
    """
    Search residues and atoms in pocket
    args: -> dictionary dataset
          -> path directiory dataset
          -> type study (surflexe or Fpocket)
    return: NONE generation file pocket residues and atoms
    """
    
    list_PDB = dico_dataset.keys ()
    for PDB_ID in list_PDB : 
        print PDB_ID, "Retrieve residues/atoms pocket"
        list_dir_pocket = pathDirectory.generateListDirPocket( PDB_ID, type_pocket, name_dataset) # path pocket by PDB directory
        path_file_complex = path_dir_dataset + PDB_ID + "/protein.pdb"
        path_filin_rsa = path_dir_dataset + PDB_ID + "/protein.rsa"
        path_filin_asa = path_dir_dataset + PDB_ID + "/protein.asa"

        if debug :
            print list_dir_pocket, "List pocket directory by PDB"
            print path_file_complex, "Path complexe"
            print path_filin_rsa, "Path file rsa"
            print path_filin_asa, "Path file asa"
        
        for path_dir_pocket in list_dir_pocket : 
            path_file_protomol_mol2 = pathDirectory.searchProtomolSurflexe (path_dir_pocket)
            if path_file_protomol_mol2 == None  : 
                if type_pocket == "surflexe": 
                    print "ERROR protmol Surflexe missing --> ", PDB_ID
                    continue
            else : 
                # convert protomol one PDB
                path_file_protomol_PDB = runOtherProg.babelConvertMol2toPDB (path_file_protomol_mol2)
            
            if type_pocket == "surflexe" : 
                # generate pocket residues and pocket atoms
                path_file_pocket_res, path_file_pocket_atom = getPocketSurflexe(path_file_protomol_PDB, path_file_complex, path_dir_pocket)
            
            elif type_pocket == "Fpocket" : 
                path_file_pocket_atom = pathDirectory.searchPocketAtomFpocket (path_dir_pocket)
                path_file_pocket_res, path_file_pocket_atom = getPocketResidues(path_file_pocket_atom, path_file_complex, path_dir_pocket, type_pocket)
            
            else :
                path_file_pocket_atom = pathDirectory.searchPocketAtomFpocket (path_dir_pocket)
                path_file_pocket_res, path_file_pocket_atom = getPocketResidues(path_file_pocket_atom, path_file_complex, path_dir_pocket, type_pocket)
                
            ####################
            #    SOLVANT ACC   #
            ####################                    
            # residues
            getAccResidues(path_filin_rsa, path_file_pocket_res)
                
            # atoms
            getAccAtom(path_filin_asa, path_file_pocket_atom, path_dir_pocket)




def getPocketResidues(path_file_pocket_atom, path_file_complex, path_dir_pocket, pocket_type_retrieve, debug = 1) : 
    """
    Retrieve pocket atoms and residues with pocket
    args: -> path file pocket atom pdb
          -> path file complex (pdb or pqr)
          -> path directory result
    return: -> NONE write files pocket_Fpocket_res.pdb and pocket_Fpocket_atom.pdb
    """
    
    path_filout_res = path_dir_pocket + "pocket_" + pocket_type_retrieve + "_res.pdb"
    path_filout_atom = path_dir_pocket + "pocket_" + pocket_type_retrieve + "_atom.pdb"
    # copy file atom for same name than surflexe
    print path_filout_atom, path_file_pocket_atom, "Pocket residues"
    os.system ("cp " + path_file_pocket_atom + " " + path_filout_atom)
    filout_res = open (path_filout_res, "w")
    
    memory_res = []
    
    pocket_parsed = PDB(path_file_pocket_atom)
    complexe_parsed = PDB(path_file_complex)
    
    for pocket_res in pocket_parsed:
        num_res_pocket = pocket_res.rNum ()
        chain_lb = pocket_res.chnLbl()
        in_list = str (num_res_pocket) + str (chain_lb)
        if debug : 
            print pocket_res.rNum (), "num pocket"
            print pocket_res.chnLbl(), "chain pocket"
            
        if in_list in memory_res : 
            continue
        else : 
            for res_complex in complexe_parsed  :
                if debug : 
                    print res_complex.rNum(), "num complexe"
                    print res_complex.chnLbl(), "chain complexe"
                    
                if res_complex.rNum() == num_res_pocket and chain_lb == res_complex.chnLbl() : 
                    filout_res.write (str(res_complex))
                    memory_res.append (in_list)
                    break
    return path_filout_res, path_filout_atom    
    

def getPocketSurflexe(path_file_protomol, path_file_PDB, path_dir_pocket, debug = 0, extend = 4):
    """
    Retrieve atom and residues (only one function -> different Stephanie script) only one loop
    args: -> path file protomol
          -> path file complexe PDB (or pqr)
          -> path directory pocket
    return: -> path file residues pocket
            -> path file atoms pocket
    """
    
    path_filout_res = path_dir_pocket + "pocket_surflexe_res.pdb"
    path_filout_atom = path_dir_pocket + "pocket_surflexe_atom.pdb"
    
#    if os.path.exists(path_filout_atom) and os.path.exists(path_filout_res) and os.path.getsize(path_filout_atom) != 0 and  os.path.getsize(path_filout_res) != 0 : 
#        return path_filout_res, path_filout_atom 

    filout_res = open(path_filout_res,'w')
    filout_atom = open (path_filout_atom, "w")
    
    memory_res = []
    memory_atom = []
    
    protomol_parsed = PDB(path_file_protomol)
    prot_parsed = PDB(path_file_PDB)
    
    for protomol_element in protomol_parsed:
        for atom_protomol in protomol_element:
            if atom_protomol.header() == "HETATM":
                xyz = atom_protomol.xyz()
                xMOL = xyz[0]
                yMOL = xyz[1]
                zMOL = xyz[2]
                for residue in prot_parsed:
                    if residue.rType() == "AMINO-ACID":
                        for atom in residue:
                            if atom.resName() != 'HOH':
                                xyzPDB = atom.xyz()
                                xPDB = xyzPDB[0]
                                yPDB = xyzPDB[1]
                                zPDB = xyzPDB[2]
                                dist = sqrt((xMOL-xPDB)**2+(yMOL-yPDB)**2+(zMOL-zPDB)**2)

                                if dist <= extend:
                                    # write residue
                                    if not int(residue.rNum()) in memory_res : 
                                        memory_res.append (int(residue.rNum()))
                                        filout_res.write(str(residue))
                                    # write atom
                                    if not int (atom.atmNum()) in memory_atom : 
                                        memory_atom.append(int(atom.atmNum()))
                                        filout_atom.write(str(atom))
    filout_res.close()
    filout_atom.close ()
    return path_filout_res, path_filout_atom


def getAccResidues( path_filin_rsa, path_filin_pocket, value_Relatif_ACC = 20.0, debug = 0 ):
    """
    Get residues exposed, exposed -> relative value
    args: -> path file rsa
          -> path file pocket
          -> value relative 
    return: NONE (write file with ACC extension)
    # rm files
    os.system ("rm " + path_filin_pocket[0:-4] + "_ACC.rsa")
    os.system ("rm " + path_filin_pocket[0:-4] + "_ACC.pdb")
    """
    if os.path.exists(path_filin_pocket[0:-4] + "_ACC.rsa") and os.path.exists(path_filin_pocket[0:-4] + "_ACC.pdb") : 
        return  

    # file out
    filout_acc_rsa = open (path_filin_pocket[0:-4] + "_ACC.rsa", "w")
    filout_acc_pocket = open (path_filin_pocket[0:-4] + "_ACC.pdb", "w")

    
    list_res_rsa_parsed = parseNACCESS.fileRSA(path_filin_rsa)
    pocket_parsed = PDB(path_filin_pocket)
    
    for pocket_res in pocket_parsed :
        res_NUM = int(pocket_res.rNum ())
        chain = pocket_res.chnLbl()
        if debug :
            print res_NUM, "Residues num"
            print chain, "Chain pocket"
        for res_rsa in list_res_rsa_parsed : 
            if debug : 
                print res_rsa["resSeq"], res_rsa["chainID"]
            if res_rsa["resSeq"] == res_NUM and res_rsa["chainID"] == chain : 
                if res_rsa["REL"] >= value_Relatif_ACC : 
                    filout_acc_pocket.write (str(pocket_res))
                    filout_acc_rsa.write (res_rsa["line"])
                    break
                    
    filout_acc_rsa.close ()
    filout_acc_pocket.close ()
    
    return path_filin_pocket[0:-4] + "_ACC.pdb"
    

def getAccAtom(path_file_asa, path_atom_pocket, value_acc = 0.0, debug = 0) :
    """
    Get residues exposed, exposed -> relative value
    args: -> path file rsa
          -> path file pocket
          -> value relative 
    return: NONE (write file with ACC extension)
    os.system ("rm " + path_atom_pocket[0:-4] + "_ACC.asa")
    os.system ("rm " + path_atom_pocket[0:-4] + "_ACC.pdb") 
    """
    
#    if os.path.exists(path_atom_pocket[0:-4] + "_ACC.asa") and os.path.exists(path_atom_pocket[0:-4] + "_ACC.pdb") : 
#        return 
    
    filout_acc_asa = open (path_atom_pocket[0:-4] + "_ACC.asa", "w")
    filout_acc_pocket = open (path_atom_pocket[0:-4] + "_ACC.pdb", "w")
    
    list_atom_rsa_parsed = parseNACCESS.fileASA(path_file_asa)
    
    
    pocket_parsed = PDB(path_atom_pocket)
    
    for pocket_res in pocket_parsed :
        for pocket_atom in pocket_res : 
            chain = pocket_atom.chnLbl()
            atom_NUM = int (pocket_atom.atmNum())
            if debug : 
                print chain, "Chain pocket"
                print atom_NUM, "atom num"
            for atom_asa in list_atom_rsa_parsed : 
                if atom_asa["atomSeq"] == atom_NUM and chain == atom_asa["chainID"] : 
                    if atom_asa["ABS"] > 0.0 :
                        if debug : 
                            print "Accessible Atom" 
                        filout_acc_pocket.write (str (pocket_atom))
                        filout_acc_asa.write (atom_asa["line"])
                        break
                    
    filout_acc_asa.close ()
    filout_acc_pocket.close ()
    
    return [path_atom_pocket[0:-4] + "_ACC.asa",path_atom_pocket[0:-4] + "_ACC.pdb"]
    
    
    

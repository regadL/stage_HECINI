"""
BORREL Alexandre
04-2012
calcul descriptor
"""
from re import search
from os import listdir, path, makedirs, system

from PDB import *
import pathDirectory
import runOtherProg
import descriptorComposition
import descriptorVolumeSize
import descriptorArea
import descriptorEnergy
import parseDescriptorRadii3
import parseFpocket
import tool


def runDescriptorGlobal (dico_dataset, path_dir_dataset, protomol_type, pocket_retrieve_type , name_dataset, compo_aa, debug = 1):
    
    list_PDB = dico_dataset.keys ()
 
    print list_PDB
    
    for PDB_ID in list_PDB : 
        print "Descriptors Protomol -> ", protomol_type, " / Pocket -> ", pocket_retrieve_type 
        list_dir_pocket = pathDirectory.generateListDirPocket(PDB_ID, pocket_retrieve_type, name_dataset)
        for path_dir_pocket in list_dir_pocket : # pocket directory
            list_files = listdir(path_dir_pocket) 
            
            # general files
            path_file_complexe = path_dir_dataset +PDB_ID +  "/protein.pdb"

            
            for file_pocket in list_files :
                if search (pocket_retrieve_type, file_pocket) : 
                    if search("_res_ACC.pdb", file_pocket) : 
                        path_pocket_acc_res = path_dir_pocket + file_pocket
                    elif search("_res.pdb", file_pocket) :
                        path_pocket_res = path_dir_pocket + file_pocket
                    elif search("atom_ACC.pdb", file_pocket) : 
                        path_pocket_acc_atom = path_dir_pocket + file_pocket
                    elif search ("atom.pdb", file_pocket): 
                        path_pocket_atom =  path_dir_pocket + file_pocket
                        
                elif search ("protomol.mol2", file_pocket) : 
                    path_surflexe_protomol_mol2 = path_dir_pocket + file_pocket
                elif search ("protomol.pdb", file_pocket) : 
                    path_surflexe_protomol_pdb = path_dir_pocket + file_pocket
                elif search("vert.pqr", file_pocket) : 
                    path_Fpocket_protomol_pdb = path_dir_pocket + file_pocket
                    
            
            if protomol_type == "Fpocket" : 
                path_Fpocket_protomol_mol2 = runOtherProg.babelPDBtoMOL2 (path_Fpocket_protomol_pdb)
                path_file_protomol_mol2 = path_Fpocket_protomol_mol2
                path_file_protomol_pdb = path_Fpocket_protomol_pdb
            elif protomol_type == "surflexe" : 
                # check if surflexe protomol exist
                if not "path_surflexe_protomol_mol2" in locals() : 
                    print "ERROR protomol missing ->", PDB_ID
                    continue 
                path_file_protomol_mol2 = path_surflexe_protomol_mol2
                path_file_protomol_pdb = path_surflexe_protomol_pdb
            else : 
                path_file_protomol_mol2 = ""
                path_file_protomol_pdb = ""
            
            if debug :
                print path_pocket_acc_res, "Pocket residues ACC"
                print path_pocket_res, "Pocket residues"
                print path_pocket_atom, "Pocket atom"
                print path_pocket_acc_atom, "Pocket acc atom"
                if protomol_type != "none" : print path_file_protomol_mol2, "Protomol (mol2)"
                if protomol_type != "none" : print path_file_protomol_pdb, "Protomol (pdb)"
                print path_file_complexe, "Complexe (pdb)"
            
            composition( path_pocket_res, path_pocket_atom, path_dir_pocket, protomol_type, pocket_retrieve_type, compo_aa = compo_aa )
            volume( path_file_protomol_mol2,path_file_protomol_pdb, path_pocket_res, path_file_complexe, path_dir_pocket, protomol_type, pocket_retrieve_type )
            surfaceArea( path_pocket_acc_atom, path_file_protomol_pdb, path_dir_pocket, path_file_complexe, PDB_ID, protomol_type, pocket_retrieve_type )
            energy( path_pocket_acc_res, path_pocket_acc_atom, path_pocket_res, path_file_complexe, path_file_protomol_mol2, path_dir_pocket, protomol_type, pocket_retrieve_type )
            radi (path_pocket_atom, path_dir_pocket, protomol_type, pocket_retrieve_type)
#             if pocket_retrieve_type != "Fpocket" :
#                 name_ligand = pathDirectory.searchLigandPDB(PDB_ID, path_dir_pocket).split ("/")[-1].split (".")[0]
#                 print name_ligand
#                 path_PDB_origin = path_dir_dataset + PDB_ID + "/" + PDB_ID + ".pdb"
#                 fpocket (name_ligand, path_PDB_origin, path_dir_pocket, protomol_type, pocket_retrieve_type)



def runDescriptor (path_pocket_atom, path_pocket_res, path_file_pocket_ACC, path_filout, path_complexe = "none", path_file_protomol = "none", name_ligand = "none"):
    
    
    filout = open (path_filout, "w")
    path_dir_temp = path.dirname(path.abspath(path_filout)) + "/temp/"
    
    ##########################################
    # Check PDB end file for RADI, END final #
    ##########################################
    tool.checkENDFinalLinePDBfile(path_pocket_atom)
    tool.checkHEADERinitialLinePDBfile(path_pocket_atom)
    
    
    ##########################
    # composition descriptor #
    ##########################
    pocketRes = PDB (path_pocket_res)
    seqRes = pocketRes.aaseq()
    
    descriptorComposition.resCompo( "pocket", seqRes, filout ) 
    descriptorComposition.resAromatic( "pocket", seqRes, filout )
    descriptorComposition.resPolaire( "pocket", seqRes, filout )
    descriptorComposition.resTiny( "pocket", seqRes, filout )
    descriptorComposition.resSmall( "pocket", seqRes, filout )
    descriptorComposition.resHydrophobe( "pocket", seqRes, filout )
    descriptorComposition.resAliphatic( "pocket", seqRes, filout )
    descriptorComposition.resNeg( "pocket", seqRes, filout )
    descriptorComposition.resPos( "pocket", seqRes, filout )
    descriptorComposition.resCharged( "pocket", seqRes, filout )
    descriptorComposition.atomes( "pocket", path_pocket_atom, filout )
    
    
    #################################
    #    volume / AREA descriptor   # -> change if need, path change in volume function descriptors
    #################################    
    
#    if path_file_protomol != "none"
#        descriptorVolumeSize.volumeChimera(filout, path_file_protomol, path_dir_result)
#        descriptorVolumeSize.compute_lambdas(path_pocket_res, path_dir_result, path_file_protomol_pdb, path_protein_pdb, filout)
    
#        descriptorArea.planarity( path_file_protomol_pdb, filout)
#        descriptorArea.narrowness( path_file_protomol_pdb, filout ) 

    #########################
    #   energy descriptor   # -> change if need, path change in volume function descriptors
    #########################  
    
    pocketRes = PDB (path_pocket_res)
    seqRes = pocketRes.aaseq()
    
    pocketAcc = PDB( path_file_pocket_ACC )
    #protein = PDB( path_file_PDB)
    seqAcc = pocketAcc.aaseq()
    descriptorComposition.charge( "pocket", seqAcc, filout)
    descriptorComposition.hydrophobicityKyte("pocket", seqRes, filout)
    
    # Function on residue files not good    
    #descriptorEnergy.ratio( path_file_pocket_res, path_result, "pocketFirst", filout )
    #descriptorEnergy.ratio( path_protomol, path_result, "protomolFirst", filout )
    path_mol2 = descriptorEnergy.ratioOriginal(path_file_pocket_ACC , "pocket", filout)
#    if not protomol_type == "none" : 
#        descriptorEnergy.ratioOriginal(path_protomol, "protomol", filout)

    ##################
    # geometry RADI #
    #################
    
    path_filout_radi = runOtherProg.runRadi (path_pocket_atom, path_dir_temp )
    path_filout_PCI = runOtherProg.runPCI (path_pocket_atom, path_dir_temp )
    
    # parse together radi and PCI descriptors -> Michel descriptor
    radi_result_parsed = parseDescriptorRadii3.fileResultByPocket (path_filout_radi, path_filout_PCI)
    
#    print radi_result_parsed
    for descriptor in radi_result_parsed.keys () :
        filout.write ("pocket_" + str (descriptor.replace(" ", "_")) + "\t" + str(radi_result_parsed[descriptor][0]) + "\n")

    system ("mv " + path_filout_radi + " " + path_dir_temp + path.basename(path_pocket_atom) + ".radi") 
    system ("mv " + path_filout_PCI + " " + path_dir_temp + path.basename(path_pocket_atom) + ".pci")  
    system ("mv " + path_mol2 + " " + path_dir_temp + path.basename(path_mol2))
    
    ###############
    #  DPOCKET    # -> no testing and move temp repertory
    ###############
    if name_ligand != "none" and path_complexe != "none" : 
        list_descriptors = ["pock_vol", "apol_as_prop", "mean_as_solv_acc", "lig_vol", "prop_polar_atm"]
        path_result = runOtherProg.runDpocket (name_ligand, path_complexe, path_dir_temp)
    
        result_parsed = parseFpocket.resultDpocket (path_result)
        print result_parsed
        for descriptor in result_parsed.keys () :
            if descriptor in list_descriptors :
                filout.write ("pocket_" + str (descriptor) + "\t" + str(result_parsed[descriptor]) + "\n")
                
                
    filout.close ()



def composition( path_pocket_res, path_pocket_atom, path_dir_pocket, protomol_type = "none", pocket_type = "none", compo_aa = 1 ):
    """
    Run composition descriptors
    args: -> path pocket residues
          -> path pocket atom
          -> path directory result
          -> protomol type
          -> pocket type
    return: NONE write in filout
    """
    print "... atomic descriptor ..."
    path_filout_descriptors =  path_dir_pocket +  "descriptor_atomic_" + protomol_type +"_" + pocket_type
    os.system ("rm " + path_filout_descriptors)
    filout = open( path_filout_descriptors, 'w' )

    pocketRes = PDB (path_pocket_res)
    seqRes = pocketRes.aaseq()


    # residues
    descriptorComposition.resGlobalCount( "pocket", seqRes, filout ) # no use in model, because volme corelated
    if compo_aa != 0 : 
        descriptorComposition.resCompo( "pocket", seqRes, filout ) 
        
    descriptorComposition.resAromatic( "pocket", seqRes, filout )
    descriptorComposition.resPolaire( "pocket", seqRes, filout )
    descriptorComposition.resTiny( "pocket", seqRes, filout )
    descriptorComposition.resHydrophobe( "pocket", seqRes, filout )
    descriptorComposition.resAliphatic( "pocket", seqRes, filout )
    descriptorComposition.resNeg( "pocket", seqRes, filout )
    descriptorComposition.resPos( "pocket", seqRes, filout )
    descriptorComposition.resCharged( "pocket", seqRes, filout )
    descriptorComposition.resPro( "pocket", seqRes, filout )
    
    # atom
    descriptorComposition.atomes( "pocket", path_pocket_atom, filout )
    
    # shpere atomisic -> ot used
    #descriptorComposition.atomistic("pocket", path_pocket_atom, filout )
    
    filout.close()           
                
                

def volume( path_file_protomol_mol2, path_file_protomol_pdb, path_file_pocket_res, path_protein_pdb, path_dir_result, protomol_type, pocket_type  ):
    
    print "... descripteurs volumiques ..."
    
    if protomol_type == "none":
        return
    path_filout = path_dir_result + "descriptor_volume_" + protomol_type + "_" + pocket_type 
    filout = open(path_filout, 'w' )
    
    descriptorVolumeSize.volumeChimera(filout, path_file_protomol_mol2, path_dir_result)
    descriptorVolumeSize.compute_lambdas(path_file_pocket_res,path_dir_result, path_file_protomol_pdb, path_protein_pdb, filout)
    
    filout.close()

def surfaceArea( path_file_pocket_acc_atom, path_file_protomol_pdb, path_result, path_complex, PDB_ID, protomol_type, pocket_type ): 
    
    print "... descripteurs surfaciques ..."
    
    path_filout = path_result + "descriptor_area_" + protomol_type + "_" + pocket_type 
    filout = open( path_filout , 'w' )
    
    descriptorArea.rugosity( path_file_pocket_acc_atom, path_complex, path_result, PDB_ID, filout )
    if not protomol_type == "none" : 
        descriptorArea.planarity( path_file_protomol_pdb, filout)
        descriptorArea.narrowness( path_file_protomol_pdb, filout ) 
    filout.close()

def energy( path_file_pocket_acc_res, path_file_pocket_acc_atom, path_file_pocket_res, path_file_PDB, path_protomol,  path_result, protomol_type, pocket_type ):
    
    print path_file_pocket_acc_atom
    print "... descripteurs energetiques ..."
    
    path_filout = path_result + "descriptor_energy_" + protomol_type + "_" + pocket_type
    os.system ("rm " + path_filout)
    filout = open(path_filout , 'w' )
    
    pocketRes = PDB (path_file_pocket_res)
    seqRes = pocketRes.aaseq()
    
    pocketAcc = PDB( path_file_pocket_acc_res )
    #protein = PDB( path_file_PDB)
    seqAcc = pocketAcc.aaseq()
    
    
    
    descriptorComposition.charge( "pocket", seqAcc, filout)
    descriptorComposition.hydrophobicityKyte("pocket", seqRes, filout)
    
    # Function on residue files not good    
    #descriptorEnergy.ratio( path_file_pocket_res, path_result, "pocketFirst", filout )
    #descriptorEnergy.ratio( path_protomol, path_result, "protomolFirst", filout )
    descriptorEnergy.ratioOriginal(path_file_pocket_acc_atom, "pocket", filout)
    if not protomol_type == "none" : 
        descriptorEnergy.ratioOriginal(path_protomol, "protomol", filout)
    
    filout.close()

def radi (path_pocket_atom, path_dir_descriptor, protomol_type, pocket_type):
    """
    Run radi and retrieve descriptors
    args:-> path atom pocket
         -> path directory descriptor
         -> protomol type
         -> pocket type
    return: NONE write file descriptor-radi
    """
    
    path_filout_radi = runOtherProg.runRadi (path_pocket_atom, path_dir_descriptor)
    path_filout_PCI = runOtherProg.runPCI (path_pocket_atom, path_dir_descriptor)
    
    # parse together radi and PCI descriptors -> Michel descriptor
    radi_result_parsed = parseDescriptorRadii3.fileResultByPocket (path_filout_radi, path_filout_PCI)
    
#    print radi_result_parsed
    
    filout_descriptor = open (path_dir_descriptor + "descriptor_radi_" + protomol_type + "_" + pocket_type, "w")
    for descriptor in radi_result_parsed.keys () :
        filout_descriptor.write ("pocket_" + str (descriptor.replace(" ", "_")) + "\t" + str(radi_result_parsed[descriptor][0]) + "\n")
    filout_descriptor.close ()
        
    
def fpocket (name_ligand, path_PDB_origin, path_dir_descriptor, protomol_type, pocket_type):  
    """
    Run dpocket for retrieve fpocket descriptors
    args: -> name ligand
          -> path PDB origin
          -> path directory descriptor
          -> protomol type
          -> pocket type
    return: -> NONE write file descriptor Fpocket
    """
    list_descriptors = ["pock_vol", "apol_as_prop", "mean_as_solv_acc", "lig_vol", "prop_polar_atm", 'Real_volume']
    filout_descriptor = open (path_dir_descriptor + "descriptor_fpocket_" + protomol_type + "_" + pocket_type, "w")
    path_result = runOtherProg.runDpocket (name_ligand, path_PDB_origin, path_dir_descriptor)
    
    result_parsed = parseFpocket.resultDpocket (path_result)
    print result_parsed
    for descriptor in result_parsed.keys () :
        if descriptor in list_descriptors :
            filout_descriptor.write ("pocket_" + str (descriptor) + "\t" + str(result_parsed[descriptor]) + "\n")
    filout_descriptor.close ()
    
    
    
    
    
    
   

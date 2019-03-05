# personal module
import pathDirectory
import downloadFile
import runOtherProg
import analysis
import writeFiles
import checkPocket
import tool
import preparePDB
import descriptor
import getResidues
import dataSet
import loadDescriptors
import identityCalcul



# Generation data set

dir_dataset = pathDirectory.dataSet ( "rugosity" )


dico_K = tool.loadFileDataset(dir_dataset + "test.txt")


### download PDB
#downloadFile.importPDB( dico_K.keys (), dir_dataset)
#
#
### manage PDB files
#preparePDB.AllPDBSepareChain(dico_K, dir_dataset)
#preparePDB.AllPDBProtonation(dico_K, dir_dataset) # -> protonation only protein.pdb
#
#### run F-pocket
#runOtherProg.globalFpocket(dir_dataset, dico_K.keys () )



##################################
##                              ##
##        Identity calcul       ##
##                              ##
##################################

# download fasta file
#downloadFile.importFasta(dico_K.keys (), dir_dataset)
#identityCalcul.pasteFastaFileGlobal(dico_K.keys (), dir_dataset)











checkPocket.selectPocketFpocket(dico_K, dir_dataset)


runOtherProg.protomolGenerationSurflexe (dico_K,  dir_dataset )

runOtherProg.globalNACCESS(dico_K, dir_dataset)







#getResidues.globalGetCompositionPocket(dico_K, dir_dataset, type_study= "Fpocket")
#
#protomol_type = "surflexe"
#pocket_retrive_type = "surflexe"
#descriptor.runDescriptor(dico_K, dir_dataset, protomol_type, pocket_retrive_type)
#
#dico_descriptors_surflexe_surflexe = loadDescriptors.generalLoadDescriptor(dico_K, protomol_type, pocket_retrive_type, write_file = 1)
###

#print "--------------------------"
#print "Calcul descriptors Fpocket"
#print "--------------------------"
#
#protomol_type = "surflexe"
#pocket_retrive_type = "Fpocket"
#descriptor.runDescriptor(dico_K, dir_dataset, protomol_type, pocket_retrive_type)
#
#dico_descriptors_fpocket_surflexe = loadDescriptors.generalLoadDescriptor(dico_K, protomol_type, pocket_retrive_type, write_file = 1)



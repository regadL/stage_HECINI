import writeFiles
import runOtherProg
import pathDirectory
import globalFonction
import loadDescriptors


def comparisonWithLigandWithoutLigand (path_directory_with_ligand, path_directory_without_ligand):
    
    
    list_files_with_ligand = pathDirectory.searchFileProba (path_directory_with_ligand)
    list_files_without_ligand = pathDirectory.searchFileProba (path_directory_without_ligand)
    
    path_directory_comparison = pathDirectory.result("comparison_estimation")
    
    print list_files_with_ligand, list_files_without_ligand
    
    makeFileForHist (list_files_with_ligand, list_files_without_ligand,  path_directory_comparison, "comparison_estimation_model", "",  debug = 0)
    
    



def makeFileForHist (list_path_files_with_ligand, list_path_files_without_ligand,  prefix_path_file, name_file, type_hist, debug = 0):
    
    
    if debug : 
        print (list_path_files_with_ligand)
        print (list_path_files_without_ligand)
    
    
    dico_score = {}
    
    parseFileScore (list_path_files_with_ligand, "with ligand", dico_score)
    parseFileScore (list_path_files_without_ligand, "without ligand", dico_score)

    path_file_hist = writeFiles.histScoreComparison (dico_score, prefix_path_file + name_file)

    if debug : 
        print dico_score

    runOtherProg.comparisonEstimationPlot (path_file_hist, type_hist)


def parseFileScore (list_path_filin, estimation_type, dico_parsed):
    
    for path_filin in list_path_filin : 
        filin = open (path_filin, "r")
        list_lines = filin.readlines ()[1:]
        filin.close ()
        
        for element_line in list_lines : 
            list_elements = element_line.strip().split (" ")
            PDB =list_elements[0].replace ("\"", "")
            value_proba = list_elements[1]
            color =  list_elements[-1]
            
            if not PDB in dico_parsed.keys () :
                dico_parsed[PDB] = {}
                dico_parsed[PDB]["color"] = color
            
            dico_parsed[PDB][estimation_type] = {}
            dico_parsed[PDB][estimation_type]["value"] =  value_proba 
        
                
    
            
    
    
def FpocketWithoutLigand  (path_directory_without_ligand, name_dataset, pocket_retrive_type, protomol_type):

    # path result
    path_directory_comparison = pathDirectory.result("comparison_with_Fpocket")

    # Fpocket SCORE
    dico_dataset = globalFonction.calculDatasetDictionary(name_dataset, 0)
    dico_Fpocket = loadDescriptors.FpocketDruggability (dico_dataset, pocket_retrive_type, name_dataset)
    
    list_files_without_ligand = pathDirectory.searchFileProba (path_directory_without_ligand)
    
    file_Fpocket = generateFileFpocketScore (dico_Fpocket, path_directory_comparison)
    
    makeFileForHist (list_files_without_ligand, file_Fpocket,  path_directory_comparison, "comparison_Fpocket", "Fpocket", debug = 1)
    
    
def generateFileFpocketScore (dico_Fpocket, path_directory) : 
    
    path_file = path_directory + "Fpocket_data"
    print dico_Fpocket
    filout = open (path_file, "w")
    for PDB in dico_Fpocket.keys () : 
        try : filout.write (dico_Fpocket[PDB]['Protein name'] + "\t" + dico_Fpocket[PDB]['Drug Score'] + "\t")
        except : continue
        if  dico_Fpocket[PDB]['druggability'] == "n" : 
            filout.write ("1\t")
        else : 
            filout.write ("2\t")
        filout.write (str (PDB) + "\n")
    filout.close ()
    return [path_file]
        
        
        
        
        
        
    
    
    
    


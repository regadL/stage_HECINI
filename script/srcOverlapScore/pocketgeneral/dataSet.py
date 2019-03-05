"""
BORREL Alexandre
04-2012
"""
import re, os
from copy import deepcopy

import parsePDB
import writeFiles
import parseWater
import runOtherProg
import globalFonction
import loadDescriptors
import PDB


list_druggable = ["Carbonic Anhydrase II", "Bifunctional dihydrofolate reductase-thymidylate synthase, DHFR-TS",
                   "Dihydrofolate Reductase","Thymidine Kinase","Sex Hormone Binding","PDE 4D",
                   "Mineralcorticoid receptor", "Glucocorticoid Receptor", "Androgen Receptor", "Progesterone Receptor",
                   "Estrogen Receptor", "Cytochrome P450 121", "CYPLI, Cytochrome P450 51, P450-14DM, P450-LIA1",
                   "HIV RT NNRT1 site", "HIV-1 Protease", "PPARg", "PPARg (other binding site)", "ADAM 33", 
                   "Acetylcholinesterase", "HMG-CoA Reductase", "DNA Gyrase B", "Plasminogen Kringle 4",
                   "P-Hydroxybenzoate Hydroxylase", "COX2, Prostaglandin G/H synthase 2", "Phospholipase A2",
                   "Beta-2-adrenergic Receptor/T4 lysozyme chimera", "Enoyl-[acyl-carrier-protein] reductase [NADH]",
                   "Coagulation factor X, Stuart Factor", "Cardiac Troponin C", "3 Alpha Hydroxysteroid Dehydrogenase",
                   "Beta-Lactamase", "PDE 5A", "THRA protein", "Cyclophilin C", "Aldose Reductase", "cAbl Kinase",
                   "CDK2", "EGFR, Epidermal growth factor receptor, ERBB", "MDM2", "P38 Map Kinase", "Renin", "HDAC",
                   "Topoisomerase II", "Xanthine oxydase", "Alpha Beta Tubulin","PTP-1B", "Carbonic Anhydrase II""Dihydrofolate Reductase",
                    "Thymidine Kinase", "Estrogen Receptor", "CYPLI, Cytochrome P450 51, P450-14DM, P450-LIA1", "HIV RT NNRT1 site",
                    "HIV-1 Protease", "Carbonic Andydrase II"]


list_no_druggable = ["cytosolic Branched Chain Aminotransferase 1, Protein ECA39",
                     "Farnesyl Diphosphate Synthase", "Glutathione S-Transferase P1-1",
                     "Glutathione S-transferase", "Nicotinate phosphoribosyltransferase",
                     "Exo-beta-D-glucosaminidase", "Phenylalanine-4-Monooxygenase",
                     "Cytochrome P450-SU1", "Cathepsin K", "Caspase 1", "Glutamine Methyltransferase, HemK",
                      "IP Phosphatase", "PH domain, Tyrosine protein kinase BTK", "Deoxyuridine nucleotide hydrolase",
                      "Glucoamylase", "Phosphatidylinositol spec. phospholipase", "Cellulase, Endoglucanase C", "Malate dehydrogenase",
                      "HIV RT Nucl. site"]

list_prodrug = ["Angiotensin converting enzyme", "human type II Inosine Monophosphate Dedhydrogenase, IMPDH", 
                "Neuraminidase", "Thrombin", "Penicillin Binding Protein", "HIV Integrase"]



l_apo128 = ["1TG3", "1ADS", "1TG9", "1R3C", "2FNK", "1M14", "2NPQ", "3DV7", "1H9N", "1CNJ", "1QE1", "1FQN", "1FQL", "1FQM", "3DVC", "3DVB",
            "1YDC", "2HDU", "2HDR", "1DLO", "2P9V", "1KE4", "1FQR", "1HEB", "1FSR", "1PW2", "1R39", "1HEA", "1F2W", "1EA5", "1L0G", "1L0F", 
            "1L0E", "1LV1", "1DFI", "1E2H", "1U13" ,"1JJB", "1CA2", "2PYC", "1L0D", "3D93", "2PKJ", "2CA2", "1FSN", "2QMV", "2B23", "2FSM", 
            "2FSL", "2FSO", "2CBC", "2CBB", "2CBA", "2CBE", "2CBD", "1H4N", "1WFC", "2H40", "1RAY", "1MUA", "1TH9", "1XGD", "2PV8", "1H5Z",
            "1PVG", "2PV5", "1HCL", "1W75", "1Q6V", "2ACR", "1CL5", "2BLS", "1RTJ", "1HPJ", "1PKR", "2OKR", "1THK", "2PVT", "1ZSC", "1R54",
            "1YO2", "1TE3", "2NWP", "1T9N", "1T9R", "3DVD", "2NXT", "1PDB", "2NXR", "2NXS", "1TB0", "1G3Z", "3DLK", "1CNI", "1CNG", "1X96", 
            "1HEC", "1HED", "2ILI", "1FB2", "2CTN", "1RZE", "1RZD", "1RZB", "1RZA", "1FSQ", "1PRG", "2J8T", "1TBT", "1RAZ", "2PTO", "2HDQ", 
            "2PTJ", "1TEU", "1TEQ", "1FR7", "2VA9", "1FR4", "1UGC", "1UGB", "1UGA", "1UGF", "1UGE", "1UGD", "4PTD", "2PTD", "1I9Y", "1EUW"]

l_apo138 = ["1TG3", "1ADS", "1TG9", "1R3C", "2FNK", "1M14", "2NPQ", "3DV7", "1H9N", "1CNJ", "1QE1", "1FQN", "1FQL", "1FQM", "3DVC", "3DVB", 
            "1YDC", "2HDU", "2HDR", "2BZ9", "1DLO", "2P9V", "1KE4", "1FQR", "1HEB", "1FSR", "1PW2", "1R39", "1HEA", "1F2W", "2R4R", "1EA5", 
            "1L0G", "1L0F", "1L0E", "2EB2", "1LV1", "1DFI", "1E2H", "1U13", "1JJB", "1CA2", "2PYC", "1L0D", "3D93", "1AJ4", "2PKJ", "2CA2", 
            "1FSN", "2QMV", "2B23", "1CAI", "2FSM", "2FSL", "2FSO", "2CBC", "2CBB", "2CBA", "2CBE", "2CBD", "1H4N", "1WFC", "2H40", "1RAY", 
            "1MUA", "1TH9", "1XGD", "1CAJ", "2PV8", "1H5Z", "1PVG", "3D92", "2PV5", "1HCL", "1W75", "1Q6V", "2ACR", "1CL5", "2BLS", "1RTJ", 
            "1HPJ", "1PKR", "2OKR", "1THK", "2PVT", "1ZSC", "1G0F", "1R54", "1YO2", "1TE3", "2NWP", "1T9N", "1T9R", "3DVD", "2NWO", "2NXT",
            "1PDB", "2NXR", "2NXS", "1TB0", "1G3Z", "3DLK", "1CNI", "1CNG", "1X96", "1HEC", "1HED", "2ILI", "2R4S", "1FB2", "2CTN", "1RZE", 
            "1RZD", "1RZB", "1RZA", "1FSQ", "1PRG", "2J8T", "1RAZ", "2PTO", "2HDQ", "2PTJ", "1TEU", "1TEQ", "1FR7", "2VA9", "1FR4", "1UGC", 
            "1UGB", "1UGA", "1UGE", "1UGD", "1SC4", "4PTD", "2PTD", "1I9Y", "1BTK", "1EU5", "1EUW"]


def formatFileDataSetKrasowki2011 ( path_file, debug=0 ):
    """
    Parse dataset file
    args: -> path file .txt
    return: -> list of PDB ID with druggability index
    """
    
    if debug : print path_file
    dic_PDB = {}
    filin = open ( path_file, "r" )
    list_lines = filin.readlines()
    filin.close()
    # regex retrieve
    regex_PDBID = re.compile ( '^[0-9][A-Za-z0-9]{3} ' )
    regex_indice_drugg = re.compile ( ' [d|n] ' )
    regex_group_data = re.compile (' [v|t] ')
    # read first part of file
    i_line = 0
    while ( list_lines [i_line] != "--------------------------------------- 4\r\n" ):
        PDB_ID = regex_PDBID.findall ( list_lines[i_line] )
        indice_drugg = regex_indice_drugg.findall ( list_lines[i_line] )
        group_data = regex_group_data.findall (list_lines[i_line])
        if debug : print indice_drugg, PDB_ID, group_data
        if PDB_ID != [] and indice_drugg !=[] and group_data !=[]:
            PDB_ID = PDB_ID[0].replace( " ", "" ).upper()
            if not PDB_ID in dic_PDB.keys () :
                dic_PDB[PDB_ID] = {}
                dic_PDB[PDB_ID]["druggability"] = indice_drugg[0].replace( " ", "" )
                dic_PDB[PDB_ID]["ligands"] = []
                if dic_PDB[PDB_ID]["druggability"] == "d" : 
                    dic_PDB[PDB_ID]["Protein name"] = list_lines [i_line].split (" d ")[0][5:]
                else : 
                    dic_PDB[PDB_ID]["Protein name"] = list_lines [i_line].split (" n ")[0][5:]
                # if not found dataset -> only training set
                if group_data != [] :
                    dic_PDB[PDB_ID]["data"] = group_data[0].replace (" ", "")
                else :
                    dic_PDB[PDB_ID]["data"] = "t"  
        i_line = i_line + 1
    return dic_PDB


def formatFileDataSchmitkeDD ( path_file, debug=0 ):
    """
    Parse file with dataset DD and del prodrug
    args: path file dataset
    return: dictionary with ligand and PDB
    """
    
    if debug : print path_file 
    dic_PDB = {}
    
    filin = open (path_file, "r")
    list_lines = filin.readlines ()
    for line_fillin in list_lines[3:] :
        line_parse = line_fillin.split (";")
        PDB_ID = line_parse [3].replace ('"',"")
        if not PDB_ID in dic_PDB.keys ()  : 
            dic_PDB  [PDB_ID] = {}
            dic_PDB  [PDB_ID] ['Type structure'] = line_parse [0].replace ('"',"")
            dic_PDB  [PDB_ID] ['Protein name'] = line_parse [1].replace ('"',"")
            dic_PDB  [PDB_ID] ['Protein function'] = line_parse [2].replace ('"',"")
            dic_PDB  [PDB_ID] ['PDB holo'] = [line_parse [4].replace ('"',"")]
            dic_PDB  [PDB_ID] ['Drug Score'] = line_parse [6].replace ('"',"")
            dic_PDB  [PDB_ID] ['Confidence'] = line_parse [7].replace ('"',"")
            dic_PDB  [PDB_ID] ['Biblio'] = line_parse [8].replace ('"',"")
            dic_PDB  [PDB_ID] ['ligands'] = [line_parse [5].replace ('"',"").upper()]
            
            dic_PDB  [PDB_ID] ['data'] = "t"
            
            if dic_PDB  [PDB_ID] ['Protein name'] in list_druggable : 
                dic_PDB  [PDB_ID] ['druggability'] = "d"
            elif dic_PDB  [PDB_ID] ['Protein name'] in list_no_druggable :
                dic_PDB  [PDB_ID] ['druggability'] = "n"
            elif dic_PDB  [PDB_ID] ['Protein name'] in list_prodrug :
                dic_PDB  [PDB_ID] ['druggability'] = "n"
                del dic_PDB  [PDB_ID]
                continue
            else : 
                print dic_PDB  [PDB_ID] ['Protein name']   
        else :
            dic_PDB  [PDB_ID] ['PDB holo'].append(line_parse [3].replace ('"',""))
            
    return dic_PDB
        
        
def compareDatasetKrasowskiSchmitke (dico_krasowski, path_file_Schmidtke, path_dir_dataset, debug = 0):
    """
    Compare dataset K and S for retrieve ligand by PDB
    debug -> write file list of PDB with ligand
    return: NONE change dico_krasowski
    """
    
    dico_schmidtke = formatFileDataSchmitkeDD(path_file_Schmidtke)
    i = 0
    list_PDB_Krasowski = dico_krasowski.keys ()
    for PDB_ID in list_PDB_Krasowski :
        path_file_pdb = path_dir_dataset + PDB_ID.upper() + "/" + PDB_ID.upper() + ".pdb"
        list_schmidtke = dico_schmidtke.keys () 
        if PDB_ID in list_schmidtke  : 
            # change list of ligand
            if debug : print dico_schmidtke[PDB_ID]['ligands']
            dico_krasowski[PDB_ID]["ligands"]= dico_schmidtke[PDB_ID]['ligands']
            i = i +1
        else : 
            list_ligands = parsePDB.retrieveListLigand(path_file_pdb) 
            dico_krasowski[PDB_ID]["ligands"] = list_ligands
            

def analysisLigandDataSet (dico_krasowski,path_dir_dataset, debug = 0):
    
    for PDB_ID in dico_krasowski.keys () :
        list_ligands = dico_krasowski[PDB_ID]["ligands"]
        path_file_PDB = path_dir_dataset + PDB_ID.upper() + "/" + PDB_ID.upper() + ".pdb"
        list_atom_parsed_PDB = parsePDB.loadCoordSectionPDB(path_file_PDB)
        if len(list_ligands) > 1 :
            for ligand in list_ligands : 
                if debug : print PDB_ID, ligand, "PDB ID and Ligand"
                list_ligand_parsed = parsePDB.retrieveLigand(list_atom_parsed_PDB, ligand)
                for list_atom_ligand in list_ligand_parsed : 
                    if ligand in dico_krasowski[PDB_ID]["ligands"] : 
                        # check ligand in protein (hetatm hooked)
                        #if parsePDB.checkLigandHooked(list_atom_parsed_PDB, list_atom_ligand) == 1:
                        #    dico_krasowski[PDB_ID]["ligands"].remove (ligand)
                        #    break
                        # check length ligand Lipinski between 20->40
                        if len (list_atom_ligand) <= 20 :
                            dico_krasowski[PDB_ID]["ligands"].remove (ligand)
                            break 
    

def manualAnalysisForKrasowskiDataSet (dico_krasowki):
    """
    MANUEL modification for select pocket !!!!!!!!
    """
    
    dico_krasowki["1U30"]["ligands"].remove ("NAG")
    dico_krasowki["2I0E"]["ligands"].remove ("SEP")
#    dico_krasowki["3F1Q"]["ligands"].remove ("ACY")
#    dico_krasowki["1V16"]["ligands"].remove ("BEN")
    dico_krasowki["1HVY"]["ligands"].remove ("BME")
    dico_krasowki["3IA4"]["ligands"].remove ("NDP")
#     dico_krasowki["2CL5"]["ligands"].remove ("BU3")
    dico_krasowki["3F1Q"]["ligands"].remove ("FMN")
    dico_krasowki["1CG0"]["ligands"].remove ("GDP")
    dico_krasowki["3ETR"]["ligands"].remove ("FAD")
    dico_krasowki["1NNC"]["ligands"].remove ("BMA")
    
    # analysis by data
    dico_krasowki["1AI2"]["ligands"].remove ("NAP")
    dico_krasowki["1AI2"]["ligands"].append ("ICA")
                                             
    dico_krasowki["1EC9"]["ligands"].remove ("IPA")
    dico_krasowki["1EC9"]["ligands"].append ("XYH")
    #dico_krasowki["1QXO"]["ligands"].remove ("EDO")
#    dico_krasowki["1X9D"]["ligands"].remove ("BU1")
    #dico_krasowki["1K8Q"]["ligands"].remove ("BOG")
    

#def selectComplex2calculDescriptor (path_file_identity): -> change protocol
#    """
#    Parsing of file identity
#    args: path file identity
#    return: file parsed in dictionary
#    """
#    
#    dico_out = {}
#    if not os.path.exists(path_file_identity) : 
#        print "ERROR -> file identity missing in dataset directory"
#        print "path file -> ", path_file_identity
#        return dico_out
#    
#    filin = open (path_file_identity, "r")
#    list_lines = filin.readlines ()
#    
#    for line_file in list_lines : 
#        list_element = line_file.strip().split (" ")
#        if not list_element[0] in dico_out : 
#            if  list_element[0][0:4]  == list_element[1][0:4] : 
#                dico_out[list_element[0]] = {}
#                dico_out[list_element[0]][list_element[1]] = list_element[2]
#        else : 
#            if  list_element[0][0:4]  == list_element[1][0:4] : 
#                dico_out[list_element[0]][list_element[1]] = list_element[2]
#    return dico_out

def formatFileDataPerola(path_file_Perola) : 
    
    dic_out = {}
    filin = open (path_file_Perola, "r")
    list_lines = filin.readlines ()
    filin.close ()
    
    for line_dataset in list_lines[1:] : 
        list_element = line_dataset.split ("\t")
        PDB_ID = list_element[3].upper ()
        ligand_ID = list_element[4].strip()
        dic_out[PDB_ID] = {}
        dic_out[PDB_ID]["ligands"] = [ligand_ID]
        dic_out[PDB_ID]["oral"] = list_element[2]
        dic_out[PDB_ID]["druggability"] = "d"
        dic_out[PDB_ID]["data"] = "t"
    return dic_out
    
    

def calculIdenticWater(PDB_ID, path_directory_dataset) : 
    """
    Calcul identity structure
    args: -> PDB ID
          -> path directory dataset
    return: -> dictionary with identity
    """
    path_dir_PDB = path_directory_dataset + PDB_ID + "/"
    list_files = os.listdir(path_dir_PDB)
    list_fasta = []
    dico_out = {}
    # retrieve fasta file path in list
    for file_directory in list_files : 
        if re.search (".fasta", file_directory) : 
            list_fasta.append (file_directory)
    
    first_key = list_fasta[0][0:-6] 
    
    
    nb_fasta = len (list_fasta)
    i = 0
    while i < nb_fasta-1 : 
        dico_out[first_key] = {}
        j = i + 1
        while j < nb_fasta : 
            if not os.path.exists( path_dir_PDB + first_key + "_" + list_fasta[j][0:-6] + ".water") :
                path_file_water = runOtherProg.water(path_dir_PDB + list_fasta[i], path_dir_PDB + list_fasta[j], path_dir_PDB + first_key + "_" + list_fasta[j][0:-6] + ".water")
            dico_out[first_key] [list_fasta[j][0:-6]] = parseWater.waterFile(path_dir_PDB + first_key + "_" + list_fasta[j][0:-6] + ".water")[3]
            j = j + 1
        i = i + 1
    
    
    
    
    
    
    
    if len (list_fasta) == 2 : 
        first_key = list_fasta[0][0:-6] 
        dico_out[first_key] = {}
        if not os.path.exists( path_dir_PDB + first_key + "_" + list_fasta[1][0:-6] + ".water") : 
            path_file_water = runOtherProg.water(path_dir_PDB + list_fasta[0], path_dir_PDB + list_fasta[1], path_dir_PDB + first_key + "_" + list_fasta[1][0:-6] + ".water")
        else : 
            path_file_water = path_dir_PDB + first_key + "_" + list_fasta[1][0:-6] + ".water"
        dico_out[first_key] [list_fasta[1][0:-6]] = parseWater.waterFile(path_file_water)[3]
    else :
        nb_fasta = len (list_fasta)
        i = 0
        while i < nb_fasta-1 : 
            first_key = list_fasta[i][0:-6]
            dico_out[first_key] = {}
            j = i + 1
            while j < nb_fasta : 
                if not os.path.exists( path_dir_PDB + first_key + "_" + list_fasta[j][0:-6] + ".water") :
                    path_file_water = runOtherProg.water(path_dir_PDB + list_fasta[i], path_dir_PDB + list_fasta[j], path_dir_PDB + first_key + "_" + list_fasta[j][0:-6] + ".water")
                dico_out[first_key] [list_fasta[j][0:-6]] = parseWater.waterFile(path_dir_PDB + first_key + "_" + list_fasta[j][0:-6] + ".water")[3]
                j = j + 1
            i = i + 1
        
    return dico_out    
    
    
def predictedByTypePocket ( file_predicted, dico_dataset, debug = 0) :
    
    if debug : print file_predicted
    dico_out = {}
    filin = open (file_predicted, "r")
    list_line = filin.readlines ()
    filin.close ()
    
    for line_file in list_line[1:] : 
        list_element = line_file.strip().split (" ")
        PDB = list_element[0].replace ("\"","")
        dico_out[PDB] = {}
        
        print PDB
        print dico_dataset[PDB]
        
        dico_out[PDB]["Protein name"] = dico_dataset[PDB]["Protein name"]
        dico_out[PDB]["value"] = str (list_element[1])
        dico_out[PDB]["drug"] = str (list_element[-1])
        
        if debug : 
            print line_file
            print PDB
            print dico_dataset[PDB]["Protein name"]
            
    writeFiles.typeProteinScore (dico_out, file_predicted)
    return file_predicted
        

def diviseTrainTest (dico_dataset, data_retrieve, debug = 0 ) : 
    
    if data_retrieve == "train" : 
        data_retrieve = "t"
    elif data_retrieve == "test" : 
        data_retrieve = "v"
    list_PDB = dico_dataset.keys ()
    nb_PDB = len (list_PDB)
    if debug : print nb_PDB, "nb IN"
    
    
    i = 0
    while i < nb_PDB : 
        if dico_dataset[list_PDB[i]]["data"] != data_retrieve:
            del dico_dataset[list_PDB[i]]
            del list_PDB[i]
            nb_PDB = nb_PDB - 1
        else : 
            i = i + 1 
    
    
    if debug : print len (list_PDB), "nb out"  
    
    

def loadDatasetWithFpocketDrugScore (name_dataset, pocket_retrieve_type):
    """
    Construct dataset dictionary with pocket drug score
    args: -> name dataset
          -> pocket retrieve type
    return: dictionary with score druggability
    """
    
    dico_dataset = globalFonction.calculDatasetDictionary(name_dataset, 0)
    dico_Fpocket = loadDescriptors.FpocketDruggability (dico_dataset, pocket_retrieve_type, name_dataset)
    
    return dico_Fpocket



def reduceDicoDataset (dico_dataset, key_reduce, value_key) :
    
    dico_out = deepcopy(dico_dataset)
    
    for PDB_ID in dico_out.keys () :
        if dico_out[PDB_ID][key_reduce] != value_key :
            print dico_out[PDB_ID] 
            del dico_out[PDB_ID]

    print dico_out
    return dico_out
    
    
    
def retrieveApoHoloFormDD (dico_Apo_Holo, l_PDB = [], debug = 0) : 
    
    dico_out = deepcopy(dico_Apo_Holo)
    
    
    for PDB_ID in dico_out.keys () : 
        if dico_out[PDB_ID]["Type structure"] != "apo structure" : 
            del dico_out[PDB_ID]
    
    
    if l_PDB != [] : 
        for PDB in dico_out.keys() : 
#             if dico_Apo_Holo[PDB][]]
            
            if not PDB in l_PDB : 
                del dico_out[PDB]
                
                print len (dico_out.keys())
                
                
    list_PDB_apo = dico_out.keys () 
    if debug : 
        print "******* APO ********"
        print list_PDB_apo
        print "********************"
    
       
    for PDB_apo in list_PDB_apo : 
        list_holo = dico_Apo_Holo[PDB_apo]['PDB holo']
        if debug : print list_holo
        if debug : print dico_Apo_Holo
        for PDB_holo in list_holo : 
            if not PDB_holo in dico_out.keys () : 
                dico_out[PDB_holo] = {}
                dico_out[PDB_holo] ['Type structure'] = "holo structure"
                dico_out[PDB_holo] ['Protein name'] = dico_Apo_Holo[PDB_apo]['Protein name'] 
                dico_out[PDB_holo] ['Protein function'] = dico_Apo_Holo[PDB_apo]['Protein function'] 
                dico_out[PDB_holo] ['PDB apo'] = PDB_apo
                dico_out[PDB_holo] ['Drug Score'] = dico_Apo_Holo[PDB_apo]['Drug Score']
                dico_out[PDB_holo] ['Confidence'] = dico_Apo_Holo[PDB_apo]['Confidence']
                dico_out[PDB_holo] ['Biblio'] = dico_Apo_Holo[PDB_apo]['Biblio'] 
                dico_out[PDB_holo] ['ligands'] = dico_Apo_Holo[PDB_apo]['ligands'] 
                try : dico_out[PDB_holo] ['druggability'] =  dico_Apo_Holo[PDB_apo]['druggability'] 
                except : 
                    if dico_out[PDB_holo] ['Drug Score'] >= 5 : 
                        dico_out[PDB_apo]['druggability'] = "d"
                        dico_out[PDB_holo] ['druggability'] = "d"
                    else : 
                        dico_out[PDB_apo]['druggability'] = "n"
                        dico_out[PDB_holo] ['druggability'] = "n"
    return dico_out
   # print len (dico_Apo_Holo.keys())
   
   
   
   
def formatFileData(path_file) : 
    
    dico_out = {}
    filin = open (path_file, "r")
    list_lines = filin.readlines ()
    filin.close ()
    
    for line_file in list_lines[1:] : 
        line_parsed = line_file.strip().split ("\t")
        # APO FORM
        dico_out[line_parsed[1]] = {}
        if line_parsed[3] == "1" : 
            dico_out[line_parsed[1]]['druggability'] =  "d"
        else : 
            dico_out[line_parsed[1]]['druggability'] =  "n"
            
        dico_out[line_parsed[1]]['ligands'] =  [line_parsed[2]] 
        dico_out[line_parsed[1]]['PDB holo'] =  [line_parsed[0]]
        dico_out[line_parsed[1]]['Type structure'] = "apo structure"
        
        # HOLO FORM
        dico_out[line_parsed[0]] = {}
        if line_parsed[3] == "1" : 
            dico_out[line_parsed[0]]['druggability'] =  "d"
        else : 
            dico_out[line_parsed[0]]['druggability'] =  "n"
        dico_out[line_parsed[0]]['ligands'] =  [line_parsed[2]] 
        dico_out[line_parsed[0]]['PDB apo'] =  [line_parsed[1]]
        dico_out[line_parsed[0]]['Type structure'] = "holo structure"
    
    return dico_out
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



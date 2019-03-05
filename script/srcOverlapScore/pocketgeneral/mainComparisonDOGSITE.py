import overlapPocket
import zipCompress
import pathDirectory
from os import listdir, path, remove, system
from re import search, compile
import runOtherProg
import main

def searchBestPocket (path_folder_result_server, path_pocket_data):
    l_file_pocket = listdir (path_folder_result_server)
    
    score_max = 0.0
    path_pocket = ""
    for path_pocket in l_file_pocket : 
        # only pocket no subpocket -> append if necessary
        if search (".pdb", path_pocket) :
            if search ("SP", path_pocket.split ("_")[-1]) :
                continue
            #print path_pocket, "PATH --------><"
            score = overlapPocket.scoreOverlap (path_folder_result_server + path_pocket, path_pocket_data)
            if score > score_max : 
                path_pocket = path_folder_result_server + path_pocket
                score_max = score

    return path_pocket, score_max



def retrieveDrugScore (path_file_score, s_PDB) :
    filin = open (path_file_score, "r")
    l_line = filin.readlines ()
    filin.close ()
    for line_file in l_line : 
        if search (s_PDB, line_file) : 
            return line_file.split (" ")[1], line_file.strip().split (" ")[-1]
    return "NA", "NA"

def retrieveScoreDogsite (path_folder, s_name_pocket) :
    l_file = listdir (path_folder)
    for  file_poc in l_file : 
        if search ("PocXls_", file_poc) or search ("SpocXls_", file_poc) and search (".txt2$", file_poc) :
            score = retrieveScoreFile (path_folder +"/" +  file_poc, s_name_pocket)
            if score != "NA" : 
                return score
    return "NA"
            

def retrieveScoreFile (path_file_pocket, s_name_pocket) : 
    filin = open (path_file_pocket, "r")
    l_lines = filin.readlines ()
    filin.close ()

    for line_file in l_lines :
        #print line_file 
        if search (s_name_pocket, line_file) :
            #print line_file.split ("\t")
            return line_file.split ("\t")[5].strip()
    return "NA"


def retrieveNamePocket (path_in):
    path_out = path_in.split ("_")[-1][:-4]
    return path_out



def manageFileDOGSITE (path_folder, PDB_ID) : 

    path_html, path_pocket, path_subpocket = searchFilePocket(path_folder, PDB_ID)
    html_parsed = parseHTMLDogsite (path_html)
    changeScore (html_parsed, path_pocket)
    changeScore (html_parsed, path_subpocket)
    

def changeScore (html_parsed, path_filin) : 
    
    filout = open (path_filin + "2", "w")
    regex = compile('[0-9]+\.[0-9]+')
    filin = open (path_filin, "r")
    list_line = filin.readlines ()
    filin.close ()
    
    for lin in list_line :
        l_float_find = regex.findall (lin)
        if len (l_float_find) < 5 : 
            continue
        volume = l_float_find[0]
        for pocket in html_parsed : 
            if volume == pocket["Volume"] : 
                split_line = lin.split ("\t")
                filout.write ("%s\t%s\t%s\t%s\t%s\t%s\n"%(split_line[0], pocket["Volume"], pocket["Surface"], pocket["Lipo surface"], pocket["Depth"], pocket["Score"]))
    filout.close ()
    

def parseHTMLDogsite (path_html) : 
    
    l_out = []
    filin = open (path_html, "r")
    content_file = filin.read ()
    filin.close ()
    regex = compile('[0-9]+\.[0-9]+')
    
    pocket_group = content_file.split ("<tr>")
    for element in pocket_group : 
        l_float_find = regex.findall (element)
        if len (l_float_find) == 5 :
            d = {}
            d["Volume"]  = l_float_find[0]
            d["Surface"]  = l_float_find[1]
            d["Lipo surface"]  = l_float_find[2]
            d["Depth"]  = l_float_find[3]
            d["Score"]  = l_float_find[4]
            l_out.append(d)
    return l_out
    
    
    
def searchFilePocket (path_folder, PDB_ID):
    
    path_html = path_folder + PDB_ID + ".html"
    if not path.isfile(path_html) :
        print path_html, "NOT DEFINE"
        return []
    else :
        l_files = listdir (path_folder + PDB_ID + "/")
        for  file_poc in l_files : 
            if search ("PocXls_", file_poc) : 
                path_pocket = path_folder + PDB_ID + "/" + file_poc
            if search ("SpocXls_", file_poc) :
                path_subpocket =  path_folder + PDB_ID + "/" + file_poc
        
    return [path_html, path_pocket, path_subpocket]
        
        

def main (path_folder, pocket_type, name_dataset, path_proba):

    # manage folder result
    zipCompress.unCompressFolder(path_folder)
    # retrieve best pocket
    l_PDB = listdir (path_folder)
    
    path_file_result = pathDirectory.result("comparisonDOGSITE/" + name_dataset + "/" + pocket_type ) + "compare.dat"
    
    filout = open (path_file_result, "w")
    
    for s_PDB in l_PDB : 
        if len (s_PDB) == 4 : # case folder with PDB
            print s_PDB
            system ("rm " + path_folder + s_PDB + "/*2")
            path_pocket_ref = pathDirectory.generatePathPocketATOM (s_PDB, name_dataset, pocket_type)
            #print path_pocket_ref
            if path_pocket_ref == "NONE" : 
                continue
            manageFileDOGSITE (path_folder, s_PDB)
            
            path_pocket_bestScore, f_overlap = searchBestPocket (path_folder + s_PDB + "/", path_pocket_ref)
                
            s_name_pocket = retrieveNamePocket(path_pocket_bestScore)
            
            drug_ref, druggable = retrieveDrugScore (path_proba, s_PDB)
            if drug_ref == "NA" : 
                continue
            drug_DOGSITE = retrieveScoreDogsite (path_folder + s_PDB, s_name_pocket)
            
            print s_PDB, "PDB----"
            print f_overlap, "overlap----"
            print drug_ref, "ref----"
            print drug_DOGSITE, "DOG----"
            
            write_line = "%s\t%s\t%s\t%s" %(s_PDB, f_overlap, drug_ref, drug_DOGSITE)
            if write_line  : 
                filout.write (str (s_PDB) + "\t" + str (f_overlap) + "\t" + str (drug_ref) + "\t" + str (drug_DOGSITE) + "\t"+ str(druggable) + "\n")
    
    filout.close ()
    
    runOtherProg.plotComparison (path_file_result)





################
##### main #####
################
path_folder = "/home/borrel/druggabilityProject/DOGSITE/"


# 
# pocket_type = "Fpocket"
# name_dataset =  "krasowski"
# path_proba = "/home/borrel/druggabilityProject/result/krasowski/Fpocket/LDA/WithoutLigand/without_ligandtest_Test_quality_predict"
# main (path_folder, pocket_type, name_dataset, path_proba)
#  
# pocket_type = "proximity"
# path_proba = "/home/borrel/druggabilityProject/result/krasowski/proximity/LDA/WithLigand/with_ligandtest_Test_quality_predict"
# main (path_folder, pocket_type, name_dataset, path_proba)

 
# name_dataset = "ApoFormClean"
# pocket_type = "Fpocket"
# path_proba = "/home/borrel/druggabilityProject/result/ApoFormClean/apo_all_pocket.data.proba"
# main (path_folder, pocket_type, name_dataset, path_proba)


#  
name_dataset = "ApoForm"
pocket_type = "Fpocket"
path_proba = "/home/borrel/druggabilityProject/result/ApoForm/apo_all_pocket.data.proba"
# main (path_folder, pocket_type, name_dataset, path_proba)


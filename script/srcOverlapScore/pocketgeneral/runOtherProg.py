"""
BORREL Alexandre
04-2012
"""

# global modules
import os, re, shutil, time

#personal module
import pathDirectory
import checkPocket
import tool
import superposeStructure
import parsePDB
import writePDBfile




fpocket = "/usr/local/bin/fpocket"
#surflex = "/home/borrel/softwares/SF-Full-Distribution-v2601/Docking/bin/surflex-dock.exe"
#chimera = "/home/borrel/softwares/chimera/bin/chimera"
naccess = "/home/lregad/Research/software/naccess/naccess"
#radi = "/home/borrel/softwares/RADI.3.0.3"
radi = "/home/lregad/Research/software/RADI.4.0.1"
PCI = "/home/lregad/Research/software/PCI"
#dpocket = "/home/borrel/softwares/fpocket2/bin/dpocket"
#soft_pdb2pqr = "/home/borrel/softwares/pdb2pqr-1.7.1a/pdb2pqr.py"
TMalign = "/home/lregad/Research/software/TMalign_src/TMalign"
MMalign = "/home/lregad/Research/software/TMalign_src/MMalign"

#cavitator = "/home/borrel/softwares/cavitator/cavitator_recompilMichel/cavitator_b"



#fpocket = "/home/alexandre/softwares/fpocket2/bin/fpocket"
#surflex = "/home/alexandre/softwares/SF-Full-Distribution-v2601/Docking/bin/surflex-dock.exe"
#chimera = "/home/alexandre/softwares/chimera/bin/chimera"
#naccess = "/home/alexandre/softwares/naccess/naccess"
#radi = "/home/alexandre/softwares/RADI.3.0.3"
#PCI = "/home/alexandre/softwares/PCI"
#dpocket = "/home/alexandre/softwares/fpocket2/bin/dpocket"
#soft_pdb2pqr = "pdb2pqr"

def globalNACCESS (dico_dataset, name_dataset):
    """
    For every protein.pdb generate asa and rsa
    args: -> dictionary with dataset
          -> path directory dataset
    return: -> NONE create NACCESS in dataset directory by PDB 
    """
    
    list_PDB = dico_dataset.keys ()
    for PDB_ID in list_PDB : 
        path_dir = pathDirectory.dataSet(name_dataset + "/" + PDB_ID)
        runNACESS (path_dir + "protein.pdb")
        

def globalFpocket ( path_dir_dataset, list_pdb, debug = 1 ):
    """
    Run Fpocket for every PDB
    args: -> path directory dataset
          -> list PDB
    return: -> for every path PDB creation directory Fpocket (out)
    """
    
    for PDB_ID in list_pdb :
        path_file_pdb = path_dir_dataset  + PDB_ID + "/protein.pdb"
        if debug : path_file_pdb
        runFpocket( path_file_pdb )


def globalCavitator ( path_dir_dataset, l_pdb, debug = 1 ):
    """
    Run Cavitato for every PDB
    args: -> path directory dataset
          -> list PDB
    return: -> for every path PDB creation directory Fpocket (out)
    """
    
    p_dir_cavitator = os.path.dirname(cavitator)
    for PDB_ID in l_pdb :
        path_file_pdb = path_dir_dataset  + PDB_ID + "/protein.pdb"
        
        # parse for cavitator -> rename residues
        l_at_parsed = parsePDB.loadCoordSectionPDB(path_file_pdb, "ATOM", debug)
        parsePDB.recountRes(l_at_parsed)
        writePDBfile.coordinateSection( p_dir_cavitator + "/101m_.pdb", l_at_parsed, recorder="ATOM", header=0)
        writePDBfile.coordinateSection( p_dir_cavitator + "/" + PDB_ID + ".pdb", l_at_parsed, recorder="ATOM", header=0)
        #shutil.copyfile(path_file_pdb, p_dir_cavitator + "/101m_.pdb") # check best name change
        
        if debug : print path_file_pdb
        runCavitator(PDB_ID, path_dir_dataset)




def pdb2pqr (path_file_pdb, exist_file = 0):
    """
    Run PDB2PQR and control exist files and check error
    args:
        -> path file pdb
    return:
        -> NONE (write file .pqr, change extension .pdb -> .pqr)
    """
    
    # check exist
    if exist_file == 1 :
        if os.path.exists(path_file_pdb[0:-4] + ".pqr") :
            print "File pqr exist !!!"
            return 1
    elif exist_file == 2 : 
        os.system ("rm " + path_file_pdb[0:-4] + ".pqr")
    
    
    # run pdb2pqr
    cmd = soft_pdb2pqr + " --with-ph=7 --chain --ffout=charmm --ff=charmm " +  path_file_pdb +" "+ path_file_pdb[0:-4] + ".pqr 2> " + path_file_pdb[0:-4] +"_pqrrun"
#    cmd = "pdb2pqr " + "--with-ph=7 --chain --ffout=charmm --ff=charmm " +  path_file_pdb +" "+ path_file_pdb[0:-4] + ".pqr 2> " + path_file_pdb[0:-4] +"_pqrrun"
    print cmd
    os.system (cmd)
    # check missing atom in protein
    if not os.path.exists(path_file_pdb[0:-4] + ".pqr") : 
        print "ERROR PDB2PQR reading error file"
        # Check error and dell if last residue is on line TER
        file_error = open (path_file_pdb[0:-4] +"_pqrrun","r")
        element_file = file_error.read ()
        if re.search ("missing backbone atoms in this protein", element_file) : 
            regex = re.compile("[0-9]+")
            res_missing = int (regex.findall (element_file) [-1])
            print res_missing
            file_PDB = open (path_file_pdb,"r")
            list_line_PDB = file_PDB.readlines ()
            file_PDB.close ()
            nb_line = len (list_line_PDB)
            i = 0
            while i < nb_line : 
                if re.search ("^ATOM", list_line_PDB[i]) : 
                    res = int (list_line_PDB[i][22:26].replace (" ", ""))
                if res == res_missing : 
                    del list_line_PDB[i]
                    nb_line = nb_line - 1
                    i = i-1
                i = i + 1
            filout = open (path_file_pdb,"w") 
            for line_PDB in list_line_PDB : 
                filout.write (line_PDB)
            filout.close ()
            
            pdb2pqr(path_file_pdb)
        else : 
            print "ERROR : pdb2pqr -> impossible generated .pqr !!!"
            return 0
    os.system ("rm " + path_file_pdb[0:-4] +".propka")
    os.system ("rm " + path_file_pdb[0:-4] +"_pqrrun")
    return 1
    

def runFpocket ( path_file_pdb, m = 3, M = 6, debug=1 ):
    
    try : os.system ("rm -r  "+ path_file_pdb[0:-4] + "_out")
    except : pass
    cmd = fpocket + " -m " + str (m) + " -M " + str (M) + " -f " + path_file_pdb
    if debug : print cmd
    os.system ( cmd )



def runCavitator (PDB_ID, p_dir_dataset):
    
    p_dir_pock = p_dir_dataset + "cavitator/"
    p_dir_cavitator = os.path.dirname(cavitator)
    
    # folder in dataset directory
    try : 
        os.makedirs( p_dir_pock, mode=0777 )
    except : 
        pass
    
    cmd = cavitator + "<" + p_dir_cavitator + "/nametarg1"
    print cmd
    os.system(cmd)
    
   
    
    try :
        os.system("rm " + p_dir_cavitator + "/*.surf")
    except : 
        pass
    
    try :
        shutil.copyfile(p_dir_cavitator + "/101m_top10.pdb", p_dir_pock + PDB_ID +"_pockets_10.pdb")
    except:
        pass
    
#     try :
    os.system("rm " + p_dir_cavitator + "/*.pdb")
#     except : 
#         pass
    

def runRscriptHisto ( filin, histo_type = "---", brk = 10, order = 0, debug=1 ):
    
    cmd = "./histograms.R " + filin + " " + histo_type + " " + str (brk) + " " + str (order)
    if debug : print cmd
    os.system ( cmd )
    

def runRscriptTtest ( path_file_descriptor1, path_file_descriptor2, debug = 0 ):
    """
    Run t-test on path descriptor 1 and path descriptor 2
    args: -> path file descriptor 1
          -> path file descriptor 2
    return: -> dictionary with values t-test
    
    NB: ne pas ecrire de fichier recupere directement la sortie standard voir module python
    """
    
    dico_ttest = {}
    dico_ttest["drugg"] = 0 
    dico_ttest["no_drugg"] = 0
    dico_ttest["p-value"] = 0
    
    cmd = "./t-testDescriptor.R " + path_file_descriptor1 + " " + path_file_descriptor2 + " 1> " + path_file_descriptor1 + "_temp"
    os.system ( cmd )
    if debug : print cmd
    filin = open ( path_file_descriptor1 + "_temp", "r" )
    list_line = filin.readlines ()
    filin.close ()
    if os.path.getsize(path_file_descriptor1 + "_temp") == 0 :
        dico_ttest["p-value"] = "NA"
        dico_ttest["drugg"] = "0"
        dico_ttest["no_drugg"] = "0"
        dico_ttest["sd_drugg"] = "0"
        dico_ttest["sd_no_drugg"] = "0"
        dico_ttest["corr"] = "NA"
    else : 
        # Pvalue -> Nan -> try
        dico_ttest["p-value"] = float(list_line[0].strip().split () [-1])
        dico_ttest["drugg"] = float(list_line[1].strip().split () [-1])
        dico_ttest["no_drugg"] = float(list_line[2].strip().split () [-1])
        dico_ttest["sd_drugg"] = float(list_line[3].strip().split () [-1])
        dico_ttest["sd_no_drugg"] = float(list_line[4].strip().split () [-1])
        
        try : 
            dico_ttest["corr"] = float(list_line[5].strip().split () [-1])
            dico_ttest["p-value_cor"] = float(list_line[6].strip().split () [-1])
        except :
            dico_ttest["corr"] = "NA"
    
    #os.system ( "rm " + path_file_descriptor1 +"_temp " + path_file_descriptor1 + " " + path_file_descriptor2 )
    if debug : print dico_ttest
    return dico_ttest


def runRACP (path_filin, name_main, debug = 1):
    """
    Run ACP with path data
    args: -> path filin
           -> factor for arrown in ACP
    return: -> run R script and draw .png
    """
    
    cmd_ACP = "./ACPDescriptor.R " + path_filin + " " + name_main + " 2> /dev/null"
    if debug : print cmd_ACP
    os.system (cmd_ACP)


def protomolGenerationSurflexe (dico_K, path_dataset, pocket_retrieve_type, name_dataset ,debug = 1 ) : 
    """
    Run for every pocket in result folder surflex on list residues
    args: -> dictionary dataset
          -> path directory dataset
    return: NONE (write file surflexe)
    NB: make log file, one day perhaps !!!!!!!!!!! 
    """
    
    for PDB_ID in dico_K.keys () : 
        print PDB_ID, "RUN surflexe"
        list_dir_pocket = pathDirectory.generateListDirPocket(PDB_ID, pocket_retrieve_type, name_dataset)
        if debug : print list_dir_pocket
        
        # transform .pqr for surflex run (just change extention)
        path_file_complexe = preparePQRforSurflexe ( path_dataset + PDB_ID + "/protein.pqr")
        for dir_pocket in list_dir_pocket :
            list_files_pocket = os.listdir (dir_pocket)
            
            for file_pocket in list_files_pocket : 
                if file_pocket[-7:] == "atm.pdb" :
                    path_file_pocket = dir_pocket + file_pocket
                    break
            
            if pocket_retrieve_type != "proximity" :     
                path_file_list_res = checkPocket.retrieveListRediduesPocket(path_file_pocket)
                # run Surflex on PDB not protonated (identic Stephanie but ?)
                surflexGenerateProtomolWithResPocket (path_file_list_res, path_file_complexe, dir_pocket)
            else : 
                path_file_ligand = pathDirectory.searchLigandPDB (PDB_ID, dir_pocket)
                surflexGenerateProtomolWithLigand(path_file_ligand, path_file_complexe, dir_pocket)


def preparePQRforSurflexe (path_complexe_pqr) : 
    """
    Change extention PQR -> PDB for surflexe run
    args: -> path file complexe pqr
    return: -> path file pdb (same file just change .pqr)
    """
    
    path_filout = path_complexe_pqr[0:-4] + "_surflexe.pdb"
    cmd_cp = "cp " + path_complexe_pqr + " " + path_filout
    os.system (cmd_cp)
    return path_filout

    
def surflexGenerateProtomolWithResPocket (path_file_list_res,path_file_pdb, path_directory_result, bloat = 2, thresh = 0.2) :
    """
    Generate protomol with Surflexe and file with residues pocket
    args: -> path file list residues
          -> path file PDB or PQR
          -> directory result
          -> option Surflexe (bloat and thresh)
    return: NONE 
    """
    list_file_result = os.listdir (path_directory_result)
    
    for file_result in list_file_result : 
        if re.search ("mol2", file_result) : # file exist but dell temp file
            # dell temp file
            os.system ("rm " + path_directory_result + "*mol2")
            os.system ("rm " + path_file_list_res[0:-4] + "-voxthresh* 2> /dev/null")
            os.system ("rm " + path_file_list_res[0:-4] + "-marked*" + " 2> /dev/null")
            os.system ("rm " + path_file_list_res[0:-4] + "-comp-*" + " 2> /dev/null")
            os.system ("rm " + path_file_list_res[0:-4] + "*probe*" + " 2> /dev/null")
            #return
   
    cmd_surflexe = surflex + " -proto_thresh " + str(thresh) + " -proto_bloat " + str(bloat) + " resproto "+path_file_list_res + " " + path_file_pdb + " " + path_file_list_res[0:-4] + " 2> /dev/null"
    print cmd_surflexe
    os.system (cmd_surflexe)
    
    # rm temp files
    os.system ("rm " + path_file_list_res[0:-4] + "-voxthresh*" + " 2> /dev/null")
    os.system ("rm " + path_file_list_res[0:-4] + "-marked*" + " 2> /dev/null")
    os.system ("rm " + path_file_list_res[0:-4] + "-comp-*" + " 2> /dev/null")
    os.system ("rm " + path_file_list_res[0:-4] + "*probe*" + " 2> /dev/null")




def surflexGenerateProtomolWithLigand (path_file_ligand, path_file_pdb, path_directory_result, bloat = 2, thresh = 0.2) :
    """
    Generate protomol with Surflexe and file with ligand (PDB)
    args: -> path file list residues
          -> path file PDB or PQR
          -> directory result
          -> option Surflexe (bloat and thresh)
    return: NONE 
    """
    
    list_file_result = os.listdir (path_directory_result)
    
    for file_result in list_file_result : 
        if re.search ("mol2", file_result) : # file exist but dell temp file
            # dell temp file
            os.system ("rm " + path_directory_result + "*mol2")
            os.system ("rm " + path_file_ligand[0:-4] + "-voxthresh* 2> /dev/null")
            os.system ("rm " + path_file_ligand[0:-4] + "-marked*" + " 2> /dev/null")
            os.system ("rm " + path_file_ligand[0:-4] + "-comp-*" + " 2> /dev/null")
            os.system ("rm " + path_file_ligand[0:-4] + "*probe*" + " 2> /dev/null")
            #return
   
    cmd_surflexe = surflex + " -proto_thresh " + str(thresh) + " -proto_bloat " + str(bloat) + " proto " + str(path_file_ligand) + " " + str(path_file_pdb) + " " + str(path_file_ligand[0:-4]) + " 2> /dev/null"
    print cmd_surflexe
    os.system (cmd_surflexe)
    
    # rm temp files
    os.system ("rm " + path_file_ligand[0:-4] + "-voxthresh*" + " 2> /dev/null")
    os.system ("rm " + path_file_ligand[0:-4] + "-marked*" + " 2> /dev/null")
    os.system ("rm " + path_file_ligand[0:-4] + "-comp-*" + " 2> /dev/null")
    os.system ("rm " + path_file_ligand[0:-4] + "*probe*" + " 2> /dev/null")



  
def runNACESS (path_file_in, probe = "", verbose = 1):
    """
    Run NACCESS with file pdb 
    args: -> filin pdb
    return: path files .asa and .rsa
    """
    
    cmd_run = naccess + " " + path_file_in
    if verbose : print cmd_run

    try :
        os.system ("rm *.asa ")
        os.system ("rm *.rsa ")
        os.system ("rm *.log ")
    except :
        pass
    
    os.system (cmd_run)
    os.system ("mv *.asa " + path_file_in[0:-4] + ".asa")
    os.system ("mv *.rsa " + path_file_in[0:-4] + ".rsa")
    os.system ("rm *.log")
    
    return [path_file_in[0:-4] + ".asa", path_file_in[0:-4] + ".rsa"]



def runNACESSWithOption (path_filin, probe_value, path_dir_result, debug = 1):
    """
    Run Naccess with option use in rugosity
    args: -> path filin
          -> value probe
          -> path directory result
    return: -> path file asa
    """
    
    print path_filin
    # path filout naccess
    name_file = path_filin.split ("/")[-1][0:-4] + str (probe_value).replace(".", "_") + ".asa"
    path_filout = path_dir_result + name_file
    
    # check if exist, only one run naccess
    if os.path.exists(path_filout) and os.path.getsize(path_filout) > 0 :
        return path_filout
     
    # move filin
    path_filin_naccess = path_filin[0:-4] + str (probe_value).replace(".", "_") + path_filin[-4:]
    cmd_cp = "cp " + path_filin + " " + path_filin_naccess
    if debug : print cmd_cp
    os.system (cmd_cp) 
    
    cmd_run_nacess = naccess + " " + path_filin_naccess + " -y -c -a -z 0.01 -p "+str(probe_value)
    if debug : print cmd_run_nacess
    os.system(cmd_run_nacess)
    os.system ("rm " + path_filin_naccess) 
    
    os.system ("mv *.asa " + path_dir_result)
    if debug : print "mv *.asa " + path_dir_result
    os.system ("rm *.log")
    
    return path_filout   
    

def  babelConvertMol2toPDB (path_file_mol2) : 
    
    path_filout = path_file_mol2[0:-4] + "pdb"
    
    if not os.path.exists(path_filout) : 
        cmd_convert = "babel -imol2 " + path_file_mol2 + " -opdb " + path_filout
        print cmd_convert
        os.system (cmd_convert)
    return path_filout
 
 
def babelPDBtoMOL2 (path_file_pdb) : 
    
    path_filout = path_file_pdb[0:-4] + ".mol2"
    cmd_convert = "babel  " + path_file_pdb + " "+ path_filout
    print cmd_convert
    os.system (cmd_convert + " 2> /dev/null")
    return path_filout
    
    
def runRcorrelation (path_file):
    """Correlation between two list in path file
    args: path file in
    return: write png
    """
    
    cmd_R = "./correlation.R " + path_file
    print cmd_R
    os.system (cmd_R)


def plotCorrelation (path_file1, path_file2, name_filout):
    
    cmd_run = "./correlation_file_general.R " + path_file1 + " " + path_file2 + " " + name_filout
    print cmd_run
    os.system (cmd_run)
    

def dotChart (path_filin):
    
    cmd_run = "./dotChartCorrelation.R " + path_filin
    print cmd_run
    os.system (cmd_run)     


def lda (path_filin1, path_filin2 = 0, path_filout=0, path_file_global = 0, leave_one_out = 0, name_barplot = '0', elimcor = "0", graph = 1, debug = 1) : 
    
    
    if leave_one_out == 1 : 
        cmd = "./lda.R " + str(path_filin1) + " 1 " + str(path_filin2)  + " " +  str(path_filout) + " " + str(name_barplot) + " " + str (elimcor) + " " + str (graph) + " " + str (path_file_global)  + " 1> " + str(path_filout)
    else : 
        cmd = "./lda.R " + path_filin1 + " 0 " + path_filin2  + " " +  path_filout + " " + name_barplot + " " + str (elimcor) + " " + str (graph) + " " + str (path_file_global)  + " 1>> " + path_filout
    
    if debug :
        print cmd
    os.system(cmd)
 
   
            
    #### problem to rm files
    #os.system ("rm " + path_filin1)
    
    #if path_filin2 != 0 : 
    #    os.system ("rm " + path_filin2)


def predictLDA (path_file_model, path_file_descriptor, path_file_result, plot = 1,  debug = 1) : 
    """prediction with model LDA"""
    
    cmd = "./predictLDA.R " + path_file_model + " " + path_file_descriptor + " " + str (plot) + " 1>> " + path_file_result 
    
    if debug == 1 : 
        print cmd
    
    os.system (cmd)
    return path_file_descriptor + ".proba", path_file_result
    
    
    
def runRhistoFpocket (prefix_name, type_analysis, begin = 2.5, end = 3.5, step = 0.5, debug = 1) : 
    
    cmd = "./FpocketParameter.R " + prefix_name + "_ " + str(begin) + " " + str(end) + " " + str(step) + " " + type_analysis
    if debug : print cmd
    os.system (cmd)
    
    
def Rcart (path_filin_global , path_filin_train, path_filin_test, elimcor = 0, debug = 1) :        
    """
    Run cart
    args: -> path file train
          -> path file test
          -> path file out
          -> value elim cor
    return: NONE png cart
    """
    
    if path_filin_test == 0 and path_filin_train == 0 : 
        cmd = "./CART.R " + str (path_filin_global) + " " + str (path_filin_train) + " " + str (path_filin_test) + " " + str(elimcor) + " 1>> " + path_filin_global + "_result"
    else : 
        cmd = "./CART.R " + str (path_filin_global) + " " + str (path_filin_train) + " " + str (path_filin_test) + " " + str(elimcor) + " 1> " + path_filin_global + "_result"
    
    if debug : print cmd
    os.system (cmd)
    
    return path_filin_global + "_result"




def RRandomForest (path_filin_global , path_filin_train, path_filin_test, elimcor = 0, debug = 1) :        
    """
    Run cart
    args: -> path file train
          -> path file test
          -> path file out
          -> value elim cor
    return: NONE png cart
    """
    
    if path_filin_test == 0 and path_filin_train == 0 : 
        cmd = "./randomForest.R " + str (path_filin_global) + " " + str (path_filin_train) + " " + str (path_filin_test) + " " + str(elimcor) + " 1>> " + path_filin_global + "_result"
    else : 
        cmd = "./randomForest.R " + str (path_filin_global) + " " + str (path_filin_train) + " " + str (path_filin_test) + " " + str(elimcor) + " 1> " + path_filin_global + "_result"
    
    if debug : print cmd
    os.system (cmd)
    
    return path_filin_global + "_result"


 
 
def Rglm (path_filin_global , path_filin_train, path_filin_test, elimcor = 0, debug = 1) :        
    """
    Run cart
    args: -> path file train
          -> path file test
          -> path file out
          -> value elim cor
    return: NONE png cart
    """
    
    if path_filin_test == 0 and path_filin_train == 0 : 
        cmd = "./Rglm.R " + str (path_filin_global) + " " + str (path_filin_train) + " " + str (path_filin_test) + " " + str(elimcor) + " 1>> " + path_filin_global + "_result"
    else : 
        cmd = "./Rglm.R " + str (path_filin_global) + " " + str (path_filin_train) + " " + str (path_filin_test) + " " + str(elimcor) + " 1> " + path_filin_global + "_result"
    
    if debug : print cmd
    os.system (cmd)
    
    return path_filin_global + "_result" 
 
 
    
    
def matrixCorr (path_filout, nb_color, debug=1) : 
    
    cmd = "./cardMatrix.R " + path_filout + " " + str (nb_color) + " 1"
    
    if debug : 
        print cmd
    os.system (cmd)


def matrixImage (path_filout, nb_color, debug=1) : 
    
    cmd = "./cardMatrix.R " + path_filout + " " + str (nb_color) + " 0"
    
    if debug : 
        print cmd
    os.system (cmd)


  

def RplotElimcor (path_filin, sensibilityCarac):
    
    cmd = "./plotLDAExtendElimcor.R " + path_filin + " " + str (sensibilityCarac) 
    print cmd
    os.system (cmd)
    
    
def MDS (path_filout) : 
    
    cmd = "./MDS.R " + path_filout
    print (cmd)
    os.system (cmd)
    
def chimeraVolume (path_script_chimera, path_dir_result) : 
    """
    Run Chimera for volume and compacity calculation
    args: -> path script chimera
          -> path directory result
    return: -> path result
    """
    
    path_filout = path_dir_result + "chimera.out"
    cmd = chimera + " --nogui < " + path_script_chimera + " >" + path_filout
    print cmd
    os.system (cmd)
    return path_filout
    

def runRadi (path_pocket_atom, path_dir_descriptor, debug = 1) : 
    
    path_file_param = writeParameterRadi (path_dir_descriptor, path_pocket_atom)
    path_filout_radi = path_dir_descriptor + "radi.out"
    
    cmd_radi = radi + " < " + path_file_param + " > " + path_filout_radi
    if debug : print cmd_radi
    os.system (cmd_radi)
    #os.system ("rm " + path_file_param)
    return path_filout_radi
    

def runPCI (path_pocket_atom, path_dir_descriptor, debug = 1):
    
    
    path_file_param = writeParameterRadi (path_dir_descriptor, path_pocket_atom)
    path_filout_PCI = path_dir_descriptor + "PCI.out"
    
    cmd_PCI = PCI + " < " + path_file_param + " > " + path_filout_PCI
    if debug : print cmd_PCI
    os.system (cmd_PCI)
    #os.system ("rm " + path_file_param)
    return path_filout_PCI


def runDpocket (name_ligand, path_PDB_origin, path_dir_descriptor, debug = 1) : 
    """
    Run Dpocket
    args: -> name ligand
          -> path PDB origin
          -> path directory
    return: file results Fpocket
    """
    
    path_file_param = writeParameterFpocket (name_ligand, path_PDB_origin, path_dir_descriptor)
    cmd = dpocket + " -f " + path_file_param + " -o " + path_dir_descriptor + "out_fpocket"
    if debug : print cmd
    os.system (cmd)
#    os.system ("rm " + path_file_param)
#    os.system ("rm " + path_dir_descriptor +  "*" + "tp.txt")
#    os.system ("rm " + path_dir_descriptor +  "*" + "tnp.txt")
    return path_dir_descriptor + "out_fpocket" + "_exp.txt"
    


def writeParameterFpocket (name_ligand, path_PDB_origin, path_dir_descriptor):
    
    path_filout = path_dir_descriptor + "parameter.fpocket"
    filout = open (path_filout, "w")
    filout.write (path_PDB_origin + " " + name_ligand + "\n")
    filout.close ()
    return path_filout 
   
       

def writeParameterRadi (path_dir_descriptor, path_pocket_atom, type_extention = "PDB") : 
    """
    Write parameter file for radi
    args: -> path directory descriptors
          -> path pocket atom
          -> type file extention
    return: -> path file parameters
    """
    
    path_filout = path_dir_descriptor + "parameter.radi"
    filout = open (path_filout, "w")
    #filout.write (type_extention + "\n" + path_pocket_atom + "\nEPSTAB C\n")
    filout.write (type_extention + "\n" + path_pocket_atom + "\nEPSTAB\n")
    filout.close ()
    return path_filout
    
    
def water(path_file_fasta1, path_file_fasta2, path_filout, gapopen = 10, gapextend = 0.5, debug = 1):
    """
    Run water from emboss with 2 fasta files and gap open option and gap extend
    args: -> file fasta 1
          -> file fasta 2
          -> gap open value
          -> gap extend value
     return: -> water file
     """
    
    cmd = "water -asequence " + path_file_fasta1 + " -bsequence " + path_file_fasta2 + " -outfile " + path_filout + " -gapopen " + str(gapopen) + " -gapextend " + str(gapextend)
    
    if debug : 
        print cmd
    os.system (cmd)   
    return path_filout


def ACPDataset (path_file_data1, path_file_data2, path_file_result, correspondance_file="0", debug = 1) : 
    
    cmd = "./ACPDataset.R " + path_file_data1 + " " + path_file_data2 + " " + correspondance_file + " " +  path_file_result
    if debug : 
        print cmd
    os.system (cmd)


def MDSDataset (path_file_data1, path_file_data2, path_file_result, correspondance_file="0", debug = 1) : 
    
    cmd = "./MDSDataset.R " + path_file_data1 + " " + path_file_data2 + " " + correspondance_file + " " +  path_file_result
    if debug : 
        print cmd
    os.system (cmd)


def ldaPrediction (path_file_train, path_file_test, path_global_train, path_global_test, path_file_correct, value_cor, path_dir, path_descriptor_signif, begin_IC = 0.35, end_IC = 0.65, debug = 1) : 
    
    cmd = "./ldamodelIC.R " + path_file_train + " " + path_file_test + " "  + path_global_train + " " + path_global_test + " " + path_file_correct + " " + path_descriptor_signif +" " +str (begin_IC) + " " + str (end_IC) + " " + str (value_cor) + " > " + path_dir + "result.txt"
    
    if debug : 
        print cmd
    os.system (cmd)
    
    
    
def SVM (path_filin1, path_filin2 = 0, path_filout = "", debug = 1):
    
    path_file_result = str(path_filout) + "_result"
    cmd = "./SVM.R " + " " + str (path_filin1) + " " + str (path_filin2) + " > " + path_file_result
    
    if debug : 
        print cmd
    
    os.system (cmd)
    return path_file_result
    
        
    
def histNameProtein (path_file_protein) : 
    
    cmd = "./histType.R " + path_file_protein
    print cmd
    os.system (cmd)



def cartPlot (path_file_descriptor, path_file_tree, debug = 0) :
    """
    Run script to draw the tree by CART
    args: -> path file descriptor
          -> path file tree
    return: NONE run R script
    """

    cmd = "./cartTree.R " + path_file_descriptor + " " + path_file_tree
    if debug : print cmd
    os.system (cmd)

def correlationDescriptor (path_filout_desc, debug = 1) : 
    """
    Draw correlation between two descriptor
    args: -> path file descriptors 1
          -> path file descriptors 2
    return: NONE draw R plot
    """
    
    cmd = "./correlationBetween2Desc.R " + path_filout_desc
    if debug : 
        print cmd
    os.system (cmd)

    
def accTrainTest (path_file_train, path_file_test, begin, step) : 
    """
    Draw plot with accuracy with train and test file
    args: -> path file train
          -> path file test
          -> begin value cor test
          -> step value cor train
    return: NONE 
    """
    
    cmd = "./parameterTestTrain.R " + path_file_train + " " + path_file_test + " " +str (begin) + " " + str (step)
    print cmd
    os.system (cmd)


def descriptorSelectionLDA (path_file_train, path_file_test, path_filout, nb_descripteur = 0, debug = 1) : 
    
    # change option select descriptor
    # test
    cmd = "./selectDescriptorLDAModel.R " + path_file_train + " " + path_file_test + " AllModel " +  str (nb_descripteur) +  " > " + path_filout   
    # correlation -> last descriptors
    #cmd = "./variableSelection.R " + path_file_train + " " + path_file_test + " " + str (cor_value) + " drugg" + " > " + path_filout   
    
    if debug : print cmd
    os.system (cmd)
    return path_filout
    
    
    
def comparisonEstimationPlot (path_file_hist, type_hist, debug = 1) :
     
    cmd = "./dotChartEstimation.R " + str(path_file_hist) + " " + str(path_file_hist)
    if debug : 
        print cmd
    
    os.system (cmd)
    

def ldaROCCurve (path_filin1, path_filin2, path_filout, debug = 1) : 
    """
    Generate ROC curve for LDA prediction -> R script
    args: -> path file trainningplotQualityMCCACC (p_filout)
           -> path file test
           -> path file png
    return: NONE create graphic
    """
    
    cmd = "./ROC_curve.R " + path_filin1 + " " + path_filin2 + " " + path_filout + "ROC.png"
    
    if debug : print cmd
    
    os.system(cmd)



def runTMalign(path_pr1, path_pr2, path_dir_out, debug = 1) : 
    
    # delet chain in PDB
    path_pr1 = tool.removeChain (path_pr1)
    path_pr2 = tool.removeChain (path_pr2)
    
    path_pr1 = superposeStructure.manageTMalign (path_pr1)
    path_pr2 = superposeStructure.manageTMalign (path_pr2)
    
    cmd_run = TMalign + " " + str (path_pr1) + " " + str (path_pr2) + " -o " + path_dir_out + "align.out -m " + path_dir_out + "matrix.out" +" > " + path_dir_out + "RMSD"
    if debug : 
        print cmd_run
    os.system(cmd_run)
    
    return [path_dir_out + "align.out", path_dir_out + "align.out_all", path_dir_out + "align.out_atm",path_dir_out + "align.out_all_atm", path_dir_out + "RMSD" ]





def runMMalign(path_pr1, path_pr2, path_dir_out, debug = 1) : 
    
    cmd_run = MMalign + " " + str (path_pr1) + " " + str (path_pr2) + " -o " + path_dir_out + "align.supp > " + path_dir_out + "RMSD"
    if debug : 
        print cmd_run
    os.system(cmd_run)
    
    return [path_dir_out + "align.out", path_dir_out + "align.supp",   path_dir_out + "RMSD" ]




   

def bartlett (path_filin , path_filout , group_variable = "drugg", debug = 1) : 
    
    cmd_run = "./bartlett.R " + path_filin + " " + group_variable + " > " + path_filout
    
    if debug : print cmd_run
    
    os.system (cmd_run)
    
    return path_filout


def descriptorSelectionSVM (path_file_train, path_file_test, path_filout, nb_cross = 10, debug = 1) : 
    
    # change option select descriptor
    # test
    cmd = "./selectDescriptorSVMModel.R " + path_file_train + " " + path_file_test + " " + str (nb_cross) + " > " + path_filout + "&"   
    # correlation
    #cmd = "./variableSelection.R " + path_file_train + " " + path_file_test + " " + str (cor_value) + " drugg" + " > " + path_filout   
    
    if debug : print cmd
#     os.system (cmd)
    return path_filout
   

def babelPDBtoSMI(path_pdb_ligand, debug = 1) : 
    
    path_filout = path_pdb_ligand[0:-4] + ".smi"
    cmd = "babel " + path_pdb_ligand + " " + path_filout
    if debug : print cmd
    os.system (cmd)
    
    return path_filout
    
    
def histogram (path_global_descriptor, path_begin, name_class, debug = 1) :
    
    cmd = "./histogram.R " + path_global_descriptor + " " + path_begin + " " + name_class  
    
    if debug : print cmd  
    
    os.system (cmd)



def histScoreByPDB (path_fillin, debug = 1) : 
    
    cmd = "./histScorePDB.R " + path_fillin
    
    if debug : 
        print cmd
        
    os.system (cmd)
    

def histoApoHoloDescRMSD(path_file_descriptor_apo, path_file_descriptor_holo, path_file_RMSD, path_directory, descriptor, debug = 1)  :
    
    
    cmd = "./histApoHoloRMSDDesc.R " + str(path_file_descriptor_apo) + " " + str(path_file_descriptor_holo) + " " + str(path_file_RMSD) +" " + str(path_directory) +   " " + str(descriptor[0])
    if debug : print cmd
    
    
    os.system (cmd)
    

def plotComparison (path_file_result, debug = 1) : 
    
    cmd = "./plotScore.R " + path_file_result + " > " + path_file_result + ".rate"
    
    if debug : print cmd
    
    os.system (cmd)
    
    
def plotScoreCriterion (path_file_result, debug = 1):
    
    cmd = "./plotCriterionScore.R " + path_file_result 
    
    if debug : print cmd
    os.system (cmd)
    
    
    
def boxplotFamily (path_proba, debug = 1) :
    
    cmd = "./boxplotFamily.R " + path_proba
    if debug : print cmd
    os.system(cmd)
    
    
def correlationPlotSO (p_filout, thresholdSO, debug = 1) : 
    
    
    cmd = "./correlationPlot.R " + p_filout + " " + str (thresholdSO)
    if debug : print cmd
    os.system (cmd)
    
    
    
def plotQualityMCCACC (p_filin, debug = 1) : 
    
    cmd = "./plotQualitySelectionModel.R " + p_filin + " 2> /dev/null"
    if debug : print cmd
    os.system (cmd)
    
    
    
def pieDescriptor (p_filin, debug = 1) : 
    
    cmd = "./pieDescriptor.R " + p_filin   
    if debug : print cmd
    os.system (cmd)



def runSaveModelLDA (path_descriptor, debug = 1):
    
    cmd = "./ldamodel.R " + path_descriptor
    if debug : print cmd
    os.system (cmd)
    
    
    
    
    



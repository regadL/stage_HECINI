"""
BORREL Alexandre
04-2012
"""
import tool
from os import system, path


def datasetWithLigand (dico_dataset, path_filout):
    """
    Write file with PDB ID, ligands and data train or test
    args: -> dictionary with dataset
          -> path filout
    return: NONE write file in path file
    """
    
    
    try : 
        filout = open (path_filout, "w")
    except : 
        print "[Error] -- open write file !"
    
    for PDB_ID in dico_dataset.keys () : 
        filout.write (str (PDB_ID) + "\t")
        if "druggability" in dico_dataset[PDB_ID].keys () :
            filout.write (dico_dataset[PDB_ID]["druggability"] + "\t")
        for ligands in dico_dataset[PDB_ID]["ligands"] : 
            filout.write (str (ligands) + " ")
        filout.write ("\t")
        filout.write (str (dico_dataset[PDB_ID]["data"])  + "\n")
    filout.close ()
        

def globalDescriptors (dico_descriptors, path_filout, type_pocket="all", debug = 0, pocket_separation = 1):
    """
    Write file descriptor for every pocket
    args: -> dictionary with descriptor calculation
          -> path filout
          -> type pocket (druggable / no druggable)
    return: NONE write file in path filout 
    """
    filout = open (path_filout, "w")
    if type_pocket == "all" : 
        list_type_pocket = dico_descriptors.keys ()
        try : list_type_pocket.remove ("data") # remove type data (training / testing)
        except : pass
        if pocket_separation == 1 : 
            if tool.checkSameOrderDescriptor (dico_descriptors) == 0 : 
                print "[ERROR -> no same descriptors]"
                return
            elif tool.checkSameOrderDescriptor (dico_descriptors) == 2 : 
                print "[ERROR -> Type pocket empty]"
                list_descriptors = tool.listDescriptor(dico_descriptors["No-Druggable"])
                print list_descriptors, "2"
                if len (list_descriptors) < 3 :
                    list_descriptors = tool.listDescriptor(dico_descriptors["Druggable"])
            else :
                list_descriptors = tool.listDescriptor(dico_descriptors["Druggable"]) 
        headerDescriptor (list_descriptors, filout)
        for type_pocket in list_type_pocket : 
            descriptorByTypePocket (dico_descriptors, type_pocket, filout)
    else : 
        list_descriptors = tool.listDescriptor(dico_descriptors[type_pocket])
        
        if debug :
            print list_descriptors, "List descriptor write file"
            print type_pocket, "Type pocket"
            print list_descriptors, "List Descriptors"
            
        headerDescriptor (list_descriptors, filout)
        descriptorByTypePocket (dico_descriptors, type_pocket, filout)
            
    filout.close ()
    return path_filout


def descriptorByTypePocket (dico_descriptors, type_pocket, filout):
    """
    Write descriptor by type pocket
    args: -> dictionary with descriptor
          -> type pocket
          -> file open write
    return: -> write in filout
    """
    
    dico_descriptors = dico_descriptors[type_pocket]
    nb_value = tool.nbValue (dico_descriptors)
    if nb_value == 0 : 
        print "[ERROR] no value in keys " + type_pocket
    for i_value in range (0,nb_value) : 
        valueDescriptor (dico_descriptors, i_value, filout, type_pocket)
    

def onePocketDescriptor (dico_desc, path_filout) : 
    filout = open (path_filout, "w")
    
    
    l_desc = dico_desc.keys ()
    headerDescriptor (l_desc, filout)
    
    for desc in l_desc : 
        filout.write (str(dico_desc[desc][0]) + "\t")
    filout.write ("0\n")
    filout.close ()





def headerDescriptor (list_descriptors, filout):
    """
    Write header descriptor files
    args: -> list descriptors
          -> file open write
    return: NONE write in filout
    """
    
    for descriptor in list_descriptors : 
        filout.write (descriptor + "\t")
    filout.write ("drugg\n") 
    


def valueDescriptor (dico_descriptors, i_value, filout, type_pocket, debug = 0 ):
    """
    Write every values by position
    args: -> dictionary with descriptors
          -> position in dictionary
          -> file open write
          -> type pocket
    return: NONE write
    """
    
    list_type_descriptor = dico_descriptors.keys ()
    list_type_descriptor.remove ("PDB")
    try : list_type_descriptor.remove ("data")
    except: pass
    filout.write (str (dico_descriptors["PDB"][i_value]) + "\t")
    for type_descriptor in list_type_descriptor : 
        for descriptor in dico_descriptors[type_descriptor].keys () : 
            filout.write (str (dico_descriptors[type_descriptor][descriptor][i_value]) + "\t")
    
    if type_pocket == "Druggable" : 
        filout.write ("1\n")
    else : 
        filout.write ("0\n")
    

def writeListValues (list_value, path_filout, NA = 0):
    """
    Write list value and append NA
    args: -> list value
          -> path file out
    return: NONE write value and close
    """
    filout = open (path_filout, "w")
    
    for value in list_value : 
        if NA == 0 and value != "NA" and value != "-nan" : 
            filout.write (str(value) + "\n")
        elif NA == 1 :
            filout.write (str(value) + "\n")
    
    filout.close ()
    

    
def specificDescriptor(dico_descriptor, descriptor_in, path_filout, debug = 1 ) : 
    """
    Write file with specific descriptor
    args: -> dictionary with every descriptor and value
          -> type descriptor or list descriptor
          -> path file out
    return: NONE write file
    """
    
    h = 0
    filout = open (path_filout, "w")
    list_type_pocket = dico_descriptor.keys ()
    list_type_pocket.remove ("data")# remove data test and train
    for type_pocket in list_type_pocket :
        nb_value = tool.nbValue (dico_descriptor[type_pocket])
        if nb_value == 0 : 
            print "[ERROR value] --> Specific"
            continue
        
        if type (descriptor_in) is str : 
            list_descriptor = dico_descriptor[type_pocket][descriptor_in].keys ()
            if h == 0 :  
                headerDescriptor (list_descriptor, filout)
                h = 1
            for i in range (0,nb_value) : 
                filout.write (dico_descriptor[type_pocket]["PDB"][i] + "\t")
                for descriptor in list_descriptor : 
                    filout.write (str (dico_descriptor[type_pocket][descriptor_in][descriptor][i]) + "\t")
                if type_pocket == "Druggable" : 
                    filout.write ("1\n")
                else : 
                    filout.write ("0\n")
        
        elif type (descriptor_in) is list : 
            list_descriptor = descriptor_in
            if h == 0 :  
                headerDescriptor (list_descriptor, filout)
                h = 1
            for i in range (0,nb_value) : 
                filout.write (dico_descriptor[type_pocket]["PDB"][i] + "\t")
                for descriptor in list_descriptor :
                    try :descriptor = descriptor.replace(" ", "_")
                    except: descriptor = descriptor
                    #if descriptor == "X._ATOM_CONVEXE" : descriptor = "%_ATOM_CONVEXE"
                    i_type = 0
                    list_type = dico_descriptor[type_pocket].keys ()
                    list_type.remove ("PDB")
                    try : 
                        list_type.remove ("data")
                    except : pass
                    if debug :
                        print descriptor
                        print dico_descriptor[type_pocket]["PDB"]
                        print type_pocket
                        
                    while not descriptor in dico_descriptor[type_pocket][list_type[i_type]].keys () :
                        i_type = i_type + 1
                        if debug : 
                            print descriptor
                            print i_type
                    filout.write (str (dico_descriptor[type_pocket][list_type[i_type]][descriptor][i]) + "\t") # case with ttest list impossible write " " in file
                if type_pocket == "Druggable" : 
                    filout.write ("1\n")
                else : 
                    filout.write ("0\n")  
    filout.close ()
    
    return path_filout
    
def specificDescriptorbyData(dico_descriptor, descriptor_in, path_filout, debug = 0 ) :  
    """
    Write values descriptor function train and test data set Krasowki
    args: -> dictionary with every descriptor
          -> type of descriptor or list descriptor
          -> path file out
    return: list type [path file train , path file test]
    NB: manage same file in train and test
    """
    
    if path.exists(path_filout + "_train") and not path.exists(path_filout + "_test"):
        return [path_filout + "_train", path_filout + "_test"]
    
    
    if debug : 
        print ("IN WRITE SPECIFIC DESCRIPTOR")
        print (dico_descriptor.keys())
        print (descriptor_in)
        print (path_filout)
        print ("____________________________")
    
    h = 0
    filout_train = open (path_filout + "_train", "w")
    filout_test = open (path_filout + "_test", "w")
    
    list_type_pocket = dico_descriptor.keys ()
    list_data = dico_descriptor["data"]
    if debug : print list_data
    
    list_type_pocket.remove ("data")
    
    for type_pocket in list_type_pocket :
        nb_value = tool.nbValue (dico_descriptor[type_pocket])
        if nb_value == 0 : 
            print "[ERROR value] --> Specific by data"
        if type (descriptor_in) is str : 
            if descriptor_in == "global" : 
                if debug : print dico_descriptor.keys ()
                list_descriptor = tool.listDescriptor(dico_descriptor["No-Druggable"])
                if h == 0 :  
                    headerDescriptor (list_descriptor, filout_train)
                    headerDescriptor (list_descriptor, filout_test)
                    h = 1
                for i in range (0,nb_value) :
                    if dico_descriptor[type_pocket]["data"][i] == "t" : 
                        filout = filout_train
                    else : 
                        filout = filout_test
                        
                    filout.write (str (dico_descriptor[type_pocket]["PDB"][i]) + "\t") 
                    for descriptor in list_descriptor :
                        list_type_descriptor = dico_descriptor[type_pocket].keys ()
                        try : list_type_descriptor.remove ("data")
                        except : pass 
                        try : list_type_descriptor.remove ("PDB")
                        except : pass
                        
                        for type_descriptor in  list_type_descriptor : 
                            #print dico_descriptor[type_pocket][type_descriptor]
                            if descriptor in dico_descriptor[type_pocket][type_descriptor].keys () : 
                                filout.write (str (dico_descriptor[type_pocket][type_descriptor][descriptor][i]) + "\t")
                                break
                    if type_pocket == "Druggable" : 
                        filout.write ("1\n")
                    else : 
                        filout.write ("0\n")
            else : 
                list_descriptor = dico_descriptor[type_pocket][descriptor_in].keys ()
                if h == 0 :  
                    headerDescriptor (list_descriptor, filout_train)
                    headerDescriptor (list_descriptor, filout_test)
                    h = 1
                
                for i in range (0,nb_value) :
                    if debug : print list_data[i]
                    if dico_descriptor[type_pocket]["data"][i]== "t" : 
                        filout = filout_train
                    else : 
                        filout = filout_test
                        
                    filout.write (str (dico_descriptor[type_pocket]["PDB"][i]) + "\t")  
                    for descriptor in list_descriptor :
                        filout.write (str (dico_descriptor[type_pocket][descriptor_in][descriptor][i]) + "\t")
                    if type_pocket == "Druggable" : 
                        filout.write ("1\n")
                    else : 
                        filout.write ("0\n")
                        
        elif type (descriptor_in) is list : 
            list_descriptor = descriptor_in
            if h == 0 :  
                headerDescriptor (list_descriptor, filout_test)
                headerDescriptor (list_descriptor, filout_train)
                h = 1
            for i in range (0,nb_value) :
                
                if debug :
                    print "-------"
                    print i, dico_descriptor[type_pocket]["data"][i]
                    print dico_descriptor[type_pocket]
                    print "-------"
                
                if dico_descriptor[type_pocket]["data"][i] == "t" : 
                    filout = filout_train
                else : 
                    filout = filout_test
                    
                filout.write (str (dico_descriptor[type_pocket]["PDB"][i]) + "\t")
                for descriptor in list_descriptor :
                    i_type = 0
                    list_type = dico_descriptor[type_pocket].keys ()
                    list_type.remove ("PDB")
                    try : list_type.remove ("data")
                    except : pass
                    if debug : 
                        print list_type
                        print dico_descriptor[type_pocket][list_type[i_type]].keys ()
                        print descriptor
                    while not descriptor in dico_descriptor[type_pocket][list_type[i_type]].keys () :
                        i_type = i_type + 1
                    filout.write (str (dico_descriptor[type_pocket][list_type[i_type]][descriptor][i]) + "\t")
                if type_pocket == "Druggable" : 
                    filout.write ("1\n")
                else : 
                    filout.write ("0\n")  
    filout_test.close ()
    filout_train.close ()
    
    # case of train data = test dataset
    if len(list(set(list_data))) == 1 : 
        system ("rm " + path_filout + "_test")
        return  [path_filout + "_train", path_filout + "_train"]
    
    
    return [path_filout + "_train", path_filout + "_test"]
     

def  colorACPFile (dico_descriptor): # change root script
        """
        Write file with color by names descriptors
        args: -> list descriptors
        return: -> path file color
        """
        l_desc_dog = ["Ellips b/a", "Ellips c/a", "Depth[A]", "Surface[A^2]", "enclosure", "Volume[A^3]", "Lipo_surf[A^2]"]
        l_desc_RADI_vol =  ["VOLUME_HULL", "SMALLEST_SIZE", "RADIUS_HULL", "DIAMETER_HULL", "SURFACE_HULL", "RADIUS_CYLINDER","C_ATOM","C_RESIDUES"] 
        l_desc_RADI_shape = ["PSI", "INERTIA_2", "INERTIA_3", "INERTIA_1", "CONVEX.SHAPE_COEFFICIENT"]
        l_desc_RADI_inter = ["PCI", "FACE", "X._ATOM_CONVEXE"]
        l_desc_Fpocket = ["Number of V. Vertices", "Mean B-factor", "Mean alpha-sphere radius", "Number of apolar alpha sphere", "Real volume (approximation)"]
        l_desc_compo = ["hydrophobic_kyte", "p_hydrophobic_residues", "p_hydrophobic_atom", "p_hyd_atom", "hydrophobicity_pocket_pocket", "p_aromatic_residues", "p_Car_atom" ,
                        "p_polar_residues", "p_aliphatic_residues", "p_Nlys_atom", "p_Ntrp_atom", "p_S_atom", "p_Otyr_atom", "p_Ooh_atom", "p_O_atom" , "p_N_atom",
                        "p_ND1_atom", "p_NE2_atom", "polarity_pocket_pocket", "p_charged_residues", "p_positive_residues", "p_negative_residues", "p_Ocoo_atom",
                        "p_Cgln_atom", "p_Ccoo_atom", "p_Carg_atom", "charge", "p_pro_residues", "p_tiny_residues",  "p_small_residues", "p_main_chain_atom", 
                        "p_side_chain_atom", "P_C_atom","p_nitrogen_atom", "p_sulfur_atom", "p_oxygen_atom",  "p_carbone_atom"]           
        
        print dico_descriptor
        path_color = "temp_color"
        file_color = open (path_color, "w")
        
        line_header = ""
        line_color = ""
#         list_type_pocket = dico_descriptor["Druggable"].keys ()
#         list_type_pocket.remove ("PDB")
#         list_type_pocket.remove ("data")
#         for type_descriptor in list_type_pocket : 
#             if type_descriptor != "fpocket" and type_descriptor != "radi" : 
#                 for descriptor in dico_descriptor["Druggable"][type_descriptor].keys () : 
#                     if descriptor == "C_ATOM" or descriptor == "C_RESIDUES" : 
#                         continue
#                     line_header = line_header + str (descriptor) + "\t"
#                     line_color = line_color + "1\t"
        
        # RADI-> VOLUME
        line_header = line_header + "\t".join(l_desc_RADI_vol) + "\t"
        line_color = line_color + "6\t6\t6\t6\t6\t6\t6\t6\t"
        
        # RADI-> SHAPE
        line_header = line_header + "\t".join(l_desc_RADI_shape) + "\t"
        line_color = line_color + "7\t7\t7\t7\t7\t"
        
        # RADI-> INTER
        line_header = line_header + "\t".join(l_desc_RADI_inter) + "\t"
        line_color = line_color + "8\t8\t8\t"
        
        # composition
        line_header = line_header + "\t".join(l_desc_compo) + "\t"
        line_color = line_color + "5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t5\t"
        
        # DogSite
        line_header = line_header + "\t".join(l_desc_dog) + "\t"
        line_color = line_color + "3\t3\t3\t3\t3\t3\t3\t"
        
        # Fpocket
        line_header = line_header + "\t".join(l_desc_Fpocket)
        line_color = line_color + "4\t4\t4\t4\t4"
        
        # END
        file_color.write (line_header + "\n")
        file_color.write (line_color + "\n")
        
        file_color.close ()
        
        return path_color   


        
def resultTtest (dico_ttest, path_filout, type_ttest, sorting = 0, option_cor = 0):
    """
    Write result means comparison
    args: -> dictionary with comparison test
          -> path filout
          -> type ttest (druggability or predicting)
          -> sorting, arrange pvalue
          -> value correlation
    return: -> NONE write in filout
    """
    
    filout = open (path_filout, "w")
    list_Pvalue = []
   
    if type_ttest == "drug" : 
        if option_cor == 1 : 
            filout.write ( "descriptors\tMdrugg\tSDdrugg\tMnodrugg\tSDnodrugg\tttestPvalue\tCoefCorr\tnb_pocket\n" )
        else : 
            filout.write ( "descriptors\tMdrugg\tSDdrugg\tMnodrugg\tSDnodrugg\tttestPvalue\tnb_pocket\n" )
    elif type_ttest == "predicted" : 
        if option_cor == 1 : 
            filout.write ( "descriptors\tMgood\tSDgodd\tMbad\tSDbad\tttestPvalue\tCoefCorr\tnb_pocket\n" )
        else : 
            filout.write ( "descriptors\tMgood\tSDgodd\tMbad\tSDbad\tttestPvalue\tnb_pocket\n" )
    
    for descriptor in dico_ttest.keys () : 
        if not sorting : 
            if option_cor == 1 : 
                filout.write ( str ( descriptor.replace( " ", "_" ) ) + "\t" + str ( dico_ttest[descriptor]["drugg"] ) + "\t" + str ( dico_ttest[descriptor]["sd_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["no_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["sd_no_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["p-value"] ) + "\t" + str (dico_ttest[descriptor]["corr"]) + "\t" + str(dico_ttest[descriptor]["nb values"]) + "\n" )
            else :
                filout.write ( str ( descriptor.replace( " ", "_" ) ) + "\t" + str ( dico_ttest[descriptor]["drugg"] ) + "\t" + str ( dico_ttest[descriptor]["sd_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["no_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["sd_no_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["p-value"] ) + "\t" + str(dico_ttest[descriptor]["nb values"]) + "\n" )
        else :
            list_Pvalue.append (dico_ttest[descriptor]["p-value"])
    
    if sorting : 
        list_Pvalue.sort()
        for pvalue in list_Pvalue :
            for descriptor in dico_ttest.keys  () : 
                if pvalue == dico_ttest[descriptor]["p-value"] : 
                    if option_cor == 1 : 
                        filout.write (str (descriptor.replace( " ", "_" )) + "\t" + str ( dico_ttest[descriptor]["drugg"] ) + "\t" + str ( dico_ttest[descriptor]["sd_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["no_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["sd_no_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["p-value"] ) + "\t" + str (dico_ttest[descriptor]["corr"]) + "\t" + str (dico_ttest[descriptor]["nb values"]) + "\n" )
                    else :
                        filout.write (str (descriptor.replace( " ", "_" ))  + "\t" + str ( dico_ttest[descriptor]["drugg"] ) + "\t" + str ( dico_ttest[descriptor]["sd_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["no_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["sd_no_drugg"] ) + "\t" + str ( dico_ttest[descriptor]["p-value"] ) + "\t" + str (dico_ttest[descriptor]["nb values"]) + "\n" )
    filout.close ()   
    
    
    
def typeProteinScore (dico_out, path_file_predicted) : 
    
    filout = open (path_file_predicted, "w")
    for PDB_ID in dico_out.keys () : 
        filout.write (str (dico_out[PDB_ID]["Protein name"]) + "\t" + str (dico_out[PDB_ID]["value"]) + "\t" +  str (dico_out[PDB_ID]["drug"]) + "\t" + str (PDB_ID) + "\n")
    filout.close ()
    return path_file_predicted    
    
    
    
def histScoreComparison (dico_score, path_filout) : 
    """
    Write file to compare two type of pocket
    args: -> dictionary with score
          -> path filout
    return: path filout (write file)
    """
    
    filout = open (path_filout, "w")
    
    for PDB in dico_score.keys () :
        filout.write (str(PDB) + "\t")
        try : 
            filout.write (str(dico_score[PDB]["with ligand"]["value"]) + "\t")
        except : 
            filout.write ("0\t")
            
        try : 
            filout.write (str(dico_score[PDB]["without ligand"]["value"]) + "\t")
        except : 
            filout.write ("0\t")
        filout.write (str(dico_score[PDB]["color"]) + "\n")
    filout.close ()
    
    return path_filout



def rateQualityFile (dico_rate, filout):
    
    
    
    filout.write ("TP -> " + str (dico_rate["TP"]) + "\n")
    filout.write ("TN -> " + str (dico_rate["TN"]) + "\n")
    filout.write ("FP -> " + str (dico_rate["FP"]) + "\n")
    filout.write ("FN -> " + str (dico_rate["FN"]) + "\n\n")
    
    filout.write ("Accuracy -> " + str (dico_rate["acc"]) + "\n")
    filout.write ("Precision -> " + str (dico_rate["pr"]) + "\n")
    filout.write ("Sensibility -> " + str (dico_rate["se"]) + "\n")
    filout.write ("Specificity -> " + str (dico_rate["sp"]) + "\n")
    filout.write ("BCR -> " + str (0.5 * (dico_rate["TP"]/(dico_rate["TP"] + dico_rate["FN"]) + dico_rate["TN"] / (dico_rate["TN"] + dico_rate["FP"]))) + "\n")


def writeCount (dico_count, path_filout) : 
    """
    Write file count with descriptors
    args : -> dictionnaty count
           -> path filout
    return: path filout
    """
    
    filout = open (path_filout, "w")
    for desc in dico_count.keys () :
        filout.write (str (desc) + "\t" + str (dico_count[desc]) + "\n")

    filout.close ()
    return path_filout


def writeCountCombi (dico_count, path_filout) : 
    """
    Write file count with descriptors
    args : -> dictionary count combinaison
           -> path filout
    return: path filout
    """    

    filout = open (path_filout, "w")
    for desc1 in dico_count.keys () :
        for desc2 in dico_count[desc1].keys () : 
            filout.write (str (desc1) + "_" + str (desc2) + "\t" + str (dico_count[desc1][desc2]) + "\n")

    filout.close ()
    return path_filout

def descImplicationLDA (dico_count, path_filout) : 
    """
    Write file count with descriptors
    args : -> dictionary count descriptor included in model
           -> path filout
    return: path filout
    """
    
    nb_des = len (dico_count.keys ()) - 1
    filout = open (path_filout, "w")
    for desc in dico_count["descriptors"] : 
        line_w = [desc]
        for i in range (0,nb_des) : 
            try : 
                line_w.append (str (dico_count[i][desc]))
            except : 
                line_w.append ("0")
                pass
        
        filout.write ("\t".join(line_w) + "\n")
        
    filout.close ()
    return path_filout
    
    
    
    
    
def corespondanceApoHolo (dico_dataset_apo_holo, path_filout) : 
    
    filout = open (path_filout, "w")
    for PDB_apo in dico_dataset_apo_holo.keys () :
        #print PDB_apo
        #print dico_dataset_apo_holo[PDB_apo].keys ()
        if "PDB holo" in dico_dataset_apo_holo[PDB_apo].keys () : 
            
            filout.write (PDB_apo + "\t" + dico_dataset_apo_holo[PDB_apo]["PDB holo"][0] + "\n")
    filout.close ()
    return path_filout
    



def listValue (dic_second, key_in, p_filout):
    
    filout = open (p_filout, "w")
    
    l_PDB = dic_second.keys ()
    
    for PDB in l_PDB : 
        for name_file in dic_second[PDB].keys () : 
            filout.write (str (dic_second[PDB][name_file][key_in]) + "\t" + str(PDB) + "\n")
    
    filout.close ()
    
    return p_filout
    
    
    
    
    
    


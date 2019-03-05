"""
BORREL Alexandre
04-2012
"""
import string
from numpy import mean

from PDB import *
import parsePDB
import tool

# Residues -> VENN Faure 2008
res_all = "ACDEFGHIKLMNPQRSTVWY"
res_pos = "HKR"
res_neg = "DE"
res_aromatic = "FYHW"
res_polar = "CDEHKNQRSTWY"
res_tiny = "ACGS"
res_hydrophobe = "GACTVLIMFWYHK"
res_aliphatic = "ILV"
res_small = "CVTGASDNP"
res_charged = "DERKH"

#res_Hacceptors = "DEHNQ"        
#res_Hdonors = "CHKNQRSTWY"



# Atoms
Atom = "CNOS"


# reference atom type (Milletti 2010) -> table 1 for aromatic carbon
dico_atom_C = {"R":["CB", "CD", "CG"], "M":["CB"], "F":["CB"], "L:":["CB","CD1", "CD2", "CG"], "W":["CB"], "D":["CB"], "K":["CB", "CE"], "H":["CB"], "V":["CB", "CG1", "CG2"], "Q":["CB", "CG", "CD"], "A":["CB"], "E":["CB"], "P":["CB", "CD", "CG"], "C":["CB"], "Y":["CB"], "N":["CB"], "I":["CB", "CD1", "CG1", "CG2"]}
dico_atom_Car = {"F":["CD1", "CD2", "CE1", "CE2", "CG", "CZ"], "W":["CD1", "CD2", "CE2", "CE3", "CG", "CH2", "CZ2", "CZ3"], "H":["CD2", "CE1"], "Y":["CD1", "CD2", "CE1", "CE2", "CG", "CZ"]}
dico_atom_Carg = {"S":["CB"], "T":["CB"]}
dico_atom_N = {"A":["N"], "C":["N"],"D":["N"],"E":["N"],"F":["N"],"G":["N"],"H":["N"],"I":["N"],"K":["N"],"L":["N"],"M":["N"],"N":["N", "ND2"],"P":["N"],"Q":["N", "NE2"],"R":["N"],"S":["N"],"T":["N"],"V":["N"],"W":["N"],"Y":["N"]}
dico_atom_ND1 = {"H":["ND1"]}
dico_atom_NE2 = {"H":["NE2"]}
dico_atom_Nlys = {"K":["NZ"]}
dico_atom_Ntrp = {"W":["NE1"]}
dico_atom_O = {"A":["O"], "C":["O"],"D":["O"],"E":["O"],"F":["O"],"G":["O"],"H":["O"],"I":["O"],"K":["O"],"L":["O"],"M":["O"],"N":["O"],"P":["O"],"Q":["O"],"R":["O"],"S":["O"],"T":["O"],"V":["O"],"W":["O"],"Y":["O"]}
dico_atom_Ocoo = {"E":["CG"]}
dico_atom_Ooh = {"S":["CB"], "T":["CB"]}
dico_atom_Otyr = {"Y":["OH"]}
dico_atom_S = {"C":["SG"]}
dico_atom_Ccoo = {"D":["CB"], "E":["CG"]}
dico_atom_Cgln = {"Q":["CG"], "N":["CB"]}
dico_atom_Hyd = {"R":["CZ"], "M":["SD"], "F":["CG", "CZ"], "L":["CG"], "W":["CE3", "CG", "CZ2"], "H":["CG"], "V":["CB"], "P":["CG"], "C":["SG"], "Y":["CG", "CZ"]}

# old version
dico_atom_hydrophobic = {"R":["CZ"], "M":["SD"], "F":["CG", "CZ"], "L":["CG"], "W":["CE3", "CG", "CZ2"], "H":["CG"], "V":["CB"], "P":["CG"], "C":["SG"], "Y":["CG", "CZ"], "I":["CB"]}
dico_Hacceptor = {"D":["OD1", "OD2"], "E":["OE1","OE2"], "H":["CG", "CD2", "NE2", "CE1", "ND1"], "N":["OD1"], "Q":["OE1"], "ASP":["OD1", "OD2"], "GLU":["OE1","OE2"], "HIS":["CG", "CD2", "NE2", "CE1", "ND1"], "ASN":["OD1"], "GLN":["OE1"]}
dico_Hdonor = {"C":["SG"], "H":["CG", "CD2", "NE2", "CE1", "ND1"], "K":["NZ"], "N":["ND2"], "Q":["NE2"], "R":["NE", "NH1", "NH2"], "S":["OG"], "T":["OG1"], "Y":["OH"], "W":["NE1"]}



# hydrophobic -> kyte doolittle
kyte_hydropathy_index = {"A": 1.8, "R":-4.5, "N":-3.5, "D":-3.5, "C":2.5, "Q":-3.5, "E":-3.5, "G":-0.4, "H":-3.2, "I":4.5, "L":3.8, "K":-3.9, "M":1.9, "F":2.8, "P":-1.6, "S":-0.8, "T":-0.7, "W":-0.9, "Y":-1.3, "V":4.2}


def resGlobalCount(flag,seq,fileout,na=0):
    """
    calcul nb_residues
    """
    
    if na==0:
        s = len(seq)
    elif na==1:
        s = "NA"
    fileout.write(flag+"_C_RESIDUES\t"+str(s)+"\n")

def resCompo(flag,seq,fileout,na = 0):
    """
    calcul proportion of amino acid
    """
    
    nb_res = float(len(seq))
    for aa in res_all:
        if na==0:
            nb = string.count(seq,aa)
        elif na==1:
            nb = "NA"
        try : fileout.write(flag+"_"+aa+"\t"+str(nb/nb_res)+"\n")
        except : fileout.write(flag+"_"+aa+"\t0\n")
        
def resAromatic(flag,seq,fileout,na = 0):
    
    nb_res = float(len(seq))
    if na==0:
        nb_aro = 0
        for aa_aro in res_aromatic:
            nb_aro += string.count(seq,aa_aro)
    elif na==1:
        nb_aro = "NA"
    try : fileout.write(flag+"_p_aromatic_residues\t"+str(nb_aro/nb_res)+"\n")
    except : fileout.write(flag+"_p_aromatic_residues\t0\n")

def resPolaire(flag,seq,fileout,na = 0):
    
    nb_res = float(len(seq))
    if na==0:
        nb_pol = 0
        for aa_pol in res_polar:
            nb_pol += string.count(seq,aa_pol)
    elif na==1:
        nb_pol = "NA"
    try : fileout.write(flag+"_p_polar_residues\t"+str(nb_pol/nb_res)+"\n")
    except : fileout.write(flag+"_p_polar_residues\t0\n")

def resTiny(flag,seq,fileout,na = 0):
    nb_res = float(len(seq))
    if na==0:
        nb_tiny = 0
        for aa_tiny in res_tiny:
            nb_tiny += string.count(seq,aa_tiny)
    elif na==1:
        nb_tiny = "NA"
    try : fileout.write(flag+"_p_tiny_residues\t"+str(nb_tiny/nb_res)+"\n")
    except : fileout.write(flag+"_p_tiny_residues\t0\n")
    


def resHydrophobe(flag,seq,fileout,na = 0):
    """
    Ratio of hydrophobic residues
    args: -> flag type of descriptor
          -> sequence
          -> filout write
    return: -> NONE (write value in filout)
    """
    
    nb_res = float(len(seq))
    if na==0:
        nb_phobe = 0
        for aa_phobe in res_hydrophobe:
            nb_phobe += string.count(seq,aa_phobe)
    elif na==1:
        nb_phobe = "NA"
    try : fileout.write(flag+"_p_hydrophobic_residues\t"+str(nb_phobe/nb_res)+"\n")
    except : fileout.write(flag+"_p_hydrophobic_residues\t0\n")

def resAliphatic(flag,seq,fileout,na= 0):
    nb_res = float(len(seq))
    if na==0:
        nb_aliph = 0
        for aa_aliph in res_aliphatic:
            nb_aliph += string.count(seq,aa_aliph)
    elif na==1:
        nb_aliph = "NA"
    try : fileout.write(flag+"_p_aliphatic_residues\t"+str(nb_aliph/nb_res)+"\n")
    except : fileout.write(flag+"_p_aliphatic_residues\t0\n")

def resNeg(flag,seq,fileout,na= 0):
    nb_res = float(len(seq))
    if na==0:
        nb_neg = 0
        for aa_neg in res_neg:
            nb_neg += string.count(seq,aa_neg)
    elif na==1:
        nb_neg = "NA"
    try : fileout.write(flag+"_p_negative_residues\t"+str(nb_neg/nb_res)+"\n")
    except : fileout.write(flag+"_p_negative_residues\t0\n")

def resPos(flag,seq,fileout,na= 0):
    nb_res = float(len(seq))
    if na==0:
        nb_pos = 0
        for aa_pos in res_pos:
            nb_pos += string.count(seq,aa_pos)
    elif na==1:
        nb_pos = "NA"
    try : fileout.write(flag+"_p_positive_residues\t"+str(nb_pos/nb_res)+"\n")
    except : fileout.write(flag+"_p_positive_residues\t0\n")

def resCharged(flag,seq,fileout,na=0):
    nb_res = float(len(seq))
    if na==0:
        nb_pos = 0
        for aa_pos in res_pos:
            nb_pos += string.count(seq,aa_pos)
        nb_neg = 0
        for aa_neg in res_neg:
            nb_neg += string.count(seq,aa_neg)
        nb_tot = nb_pos+nb_neg
    elif na==1:
        nb_tot = "NA"
    try : fileout.write(flag+"_p_charged_residues\t"+str(nb_tot/nb_res)+"\n")
    except : fileout.write(flag+"_p_charged_residues\t0\n")


def resSmall(flag,seq,fileout,na=0):
    
    nb_res = float(len(seq))
    if na==0:
        nb_small = 0
        for aa_small in res_small:
            nb_small += string.count(seq,aa_small)
    elif na==1:
        nb_small = "NA"
    try : fileout.write(flag+"_p_small_residues\t"+str(nb_small/nb_res)+"\n")
    except : fileout.write(flag+"_p_small_residues\t0\n")


def resPro(flag,seq,fileout,na=0):
    
    nb_res = float(len(seq))
    nb_pro = 0
    if na==0:
        nb_pro += string.count(seq,"P")
    elif na==1:
        nb_pro = "NA"
        
    try : fileout.write(flag+"_p_pro_residues\t"+str(nb_pro/nb_res)+"\n")
    except : fileout.write(flag+"_p_pro_residues\t0\n")


def atomes(flag, filin, filout):
    """
    Count atome
    args: -> flag filout
          -> filin, file PDB
          -> filout, file write descriptors
    return: NONE
    
    """
    
    nb_atom = 0.0
    nb_H = 0
    nb_C = 0
    nb_N = 0
    nb_O = 0
    nb_S = 0
    nb_main_chain = 0
    
    pocket_parsed = parsePDB.loadCoordSectionPDB(filin, "ATOM")
    
    C = countAtom(dico_atom_C, pocket_parsed)
    Car = countAtom(dico_atom_Car, pocket_parsed)
    Carg = countAtom(dico_atom_Carg, pocket_parsed)
    N = countAtom(dico_atom_N, pocket_parsed)
    ND1 = countAtom(dico_atom_ND1, pocket_parsed)
    NE2 = countAtom(dico_atom_NE2, pocket_parsed)
    Nlys = countAtom(dico_atom_Nlys, pocket_parsed)
    Ntrp = countAtom(dico_atom_Ntrp, pocket_parsed)
    O = countAtom(dico_atom_O, pocket_parsed)
    Ocoo = countAtom(dico_atom_Ocoo, pocket_parsed)
    Ooh = countAtom(dico_atom_Ooh, pocket_parsed)
    Otyr = countAtom(dico_atom_Otyr, pocket_parsed)
    S = countAtom(dico_atom_S, pocket_parsed)
    Ccoo = countAtom(dico_atom_Ccoo, pocket_parsed)
    Cgln = countAtom(dico_atom_Cgln, pocket_parsed)
    Hyd = countAtom(dico_atom_Hyd, pocket_parsed)
    
    Hacceptor = countAtom(dico_Hacceptor, pocket_parsed)
    Hdonor = countAtom(dico_Hdonor, pocket_parsed)
    hydrophobic = countAtom(dico_atom_hydrophobic, pocket_parsed)
    
    for atom in pocket_parsed :
        nb_atom = nb_atom + 1
        element = atom["element"]
        name_atom = atom["name"]
        if element == "C" : 
            nb_C = nb_C + 1
        elif element == "N" : 
            nb_N = nb_N + 1
        elif element  == "S" : 
            nb_S = nb_S + 1
        elif element == "O" : 
            nb_O  = nb_O + 1
        elif element == "H" : 
            nb_H = nb_H + 1
            
        if name_atom == "CA" or name_atom == "O" or name_atom == "C" or name_atom == "N" or name_atom == "H" or name_atom == "HA" :
            nb_main_chain = nb_main_chain + 1
            
    filout.write (flag + "_C_ATOM\t" + str (nb_atom) + "\n")
    filout.write (flag + "_p_side_chain_atom\t" + str (1-(nb_main_chain / nb_atom)) + "\n")
    filout.write (flag + "_p_main_chain_atom\t" + str (nb_main_chain / nb_atom) + "\n")
    filout.write (flag + "_p_sulfur_atom\t" + str (nb_S / nb_atom) + "\n")
    filout.write (flag + "_p_carbone_atom\t" + str (nb_C / nb_atom) + "\n")
    filout.write (flag + "_p_nitrogen_atom\t" + str (nb_N / nb_atom) + "\n")
    filout.write (flag + "_p_oxygen_atom\t" + str (nb_O / nb_atom) + "\n")

    # old atomistic
    #filout.write (flag + "_p_hbond_acceptor_atom\t" + str (Hacceptor / nb_atom) + "\n")
    #filout.write (flag + "_p_hbond_donor_atom\t" + str (Hdonor / nb_atom) + "\n")
    filout.write (flag + "_p_hydrophobic_atom\t" + str (hydrophobic / nb_atom) + "\n")
    
    # MILETTI
    filout.write (flag + "_P_C_atom\t" + str (C / nb_atom) + "\n")
    filout.write (flag + "_p_Car_atom\t" + str (Car / nb_atom) + "\n")
    filout.write (flag + "_p_Carg_atom\t" + str (Carg / nb_atom) + "\n")
    filout.write (flag + "_p_N_atom\t" + str (N / nb_atom) + "\n")
    filout.write (flag + "_p_ND1_atom\t" + str (ND1 / nb_atom) + "\n")
    filout.write (flag + "_p_NE2_atom\t" + str (NE2 / nb_atom) + "\n")
    filout.write (flag + "_p_Nlys_atom\t" + str (Nlys / nb_atom) + "\n")
    filout.write (flag + "_p_Ntrp_atom\t" + str (Ntrp / nb_atom) + "\n")
    filout.write (flag + "_p_O_atom\t" + str (O / nb_atom) + "\n")
    filout.write (flag + "_p_Ocoo_atom\t" + str (Ocoo / nb_atom) + "\n")
    filout.write (flag + "_p_Ooh_atom\t" + str (Ooh / nb_atom) + "\n")
    filout.write (flag + "_p_Otyr_atom\t" + str (Otyr / nb_atom) + "\n")
    filout.write (flag + "_p_S_atom\t" + str (S / nb_atom) + "\n")
    filout.write (flag + "_p_Ccoo_atom\t" + str (Ccoo / nb_atom) + "\n")
    filout.write (flag + "_p_Cgln_atom\t" + str (Cgln / nb_atom) + "\n")
    filout.write (flag + "_p_hyd_atom\t" + str (Hyd / nb_atom) + "\n")
    
    


def atomisticSphere (flag, filin, filout, max_distance = 15, analysis = 1, atom_central = "mean_point",  debug = 1):
    """Distance calcul des cercle a partir du point barycentre calcule par RADI mais non conserve car pb pas bcp d'info"""
 
    list_atom_pocket = parsePDB.loadCoordSectionPDB(filin)
    dico_stock_count = tool.generateStructCompositionAtomistic (max_distance, 3)
    
    if atom_central == "mean_point" : 
        central_point = generateMeansPointPocket (list_atom_pocket)
    # else append barycenter pocket calculated by RADI
    
    for atom in list_atom_pocket : 
        distance = parsePDB.distanceTwoatoms(central_point, atom)
       # print distance
        element = atom["element"]
        name_atom = atom["name"]
        residue = tool.transformAA(atom["resName"])
    
        for distance_key in dico_stock_count.keys() : 
            if distance <= distance_key or distance > max_distance : 
                dico_stock_count  [distance_key] ["atom"] = dico_stock_count  [distance_key] ["atom"] + 1
                if element == "C" : 
                    dico_stock_count  [distance_key] ["carbon"] = dico_stock_count  [distance_key] ["carbon"] + 1
                elif element == "N" : 
                    dico_stock_count  [distance_key] ["nitrogen"] = dico_stock_count  [distance_key] ["nitrogen"] + 1
                elif element  == "S" : 
                    dico_stock_count  [distance_key] ["sulfur"] = dico_stock_count  [distance_key] ["sulfur"] + 1
                elif element == "O" : 
                    dico_stock_count  [distance_key] ["oxygen"] = dico_stock_count  [distance_key] ["oxygen"] + 1
                elif element == "H" : 
                    dico_stock_count  [distance_key] ["hydrogen"] = dico_stock_count  [distance_key] ["hydrogen"] + 1
                
                if residue in dico_Hacceptor.keys () : 
                    if name_atom in dico_Hacceptor[residue] : 
                        dico_stock_count  [distance_key] ["hbond_acceptor"] = dico_stock_count  [distance_key] ["hbond_acceptor"] + 1
                
                if residue in dico_atom_Car : 
                    if name_atom in dico_atom_Car[residue] : 
                        dico_stock_count  [distance_key] ["aromatic"] = dico_stock_count  [distance_key] ["aromatic"] + 1
                
                if residue in dico_atom_hydrophobic : 
                    if name_atom in dico_atom_hydrophobic[residue] : 
                        dico_stock_count  [distance_key] ["hydrophobic"] = dico_stock_count  [distance_key] ["hydrophobic"] + 1
                        
                if residue in dico_atom_Carg : 
                    if name_atom in dico_atom_Carg[residue] : 
                        dico_stock_count  [distance_key] ["alcool"] = dico_stock_count  [distance_key] ["alcool"] + 1
                        
                
                if residue in dico_Hdonor.keys () : 
                    if name_atom in dico_Hdonor[residue] : 
                        dico_stock_count  [distance_key] ["hbond_donor"] = dico_stock_count  [distance_key] ["hbond_donor"] + 1
    
                if name_atom == "CA" or name_atom == "O" or name_atom == "C" or name_atom == "N" or name_atom == "H" or name_atom == "HA" :
                    dico_stock_count  [distance_key] ["main_chain"] = dico_stock_count  [distance_key] ["main_chain"] + 1
                else : 
                    dico_stock_count  [distance_key] ["side_chain"] = dico_stock_count  [distance_key] ["side_chain"] + 1
    
    for distance_key in dico_stock_count.keys () : 
        nb_atom = float(dico_stock_count  [distance_key] ["atom"])
        if nb_atom == 0 : 
            filout.write (flag + "_atom_" + str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_side_chain_"+ str(distance_key) + "\t" + "0" + "\n")
            filout.write (flag + "_main_chain_" + str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_sulfur_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_carbone_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_nitrogen_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_oxygen_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_hydrogen_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_hbond_acceptor_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_hbond_donor_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_alcool_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_hydrophobic_"+ str(distance_key) +"\t" + "0" + "\n")
            filout.write (flag + "_aromatic_"+ str(distance_key) +"\t" + "0" + "\n")
            
        else : 
            filout.write (flag + "_atom_" + str(distance_key) +"\t" + str(nb_atom) + "\n")
            filout.write (flag + "_side_chain_"+ str(distance_key) + "\t" + str (dico_stock_count  [distance_key] ["side_chain"] / nb_atom) + "\n")
            filout.write (flag + "_main_chain_" + str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["main_chain"] / nb_atom) + "\n")
            filout.write (flag + "_sulfur_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["sulfur"] / nb_atom) + "\n")
            filout.write (flag + "_carbone_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["carbon"] / nb_atom) + "\n")
            filout.write (flag + "_nitrogen_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["nitrogen"] / nb_atom) + "\n")
            filout.write (flag + "_oxygen_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["oxygen"] / nb_atom) + "\n")
            filout.write (flag + "_hydrogen_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["hydrogen"] / nb_atom) + "\n")
            filout.write (flag + "_hbond_acceptor_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["hbond_acceptor"] / nb_atom) + "\n")
            filout.write (flag + "_hbond_donor_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["hbond_donor"] / nb_atom) + "\n")
            filout.write (flag + "_alcool_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["alcool"] / nb_atom) + "\n")
            filout.write (flag + "_hydrophobic_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["hydrophobic"] / nb_atom) + "\n")
            filout.write (flag + "_aromatic_"+ str(distance_key) +"\t" + str (dico_stock_count  [distance_key] ["aromatic"] / nb_atom) + "\n")
    
def generateMeansPointPocket (list_atom_pocket) : 
    
    atom_out = {}
    list_x = []
    list_y = []
    list_z = []
    for atom in list_atom_pocket : 
        list_x.append (atom["x"])
        list_y.append (atom["y"])
        list_z.append (atom["z"])
        
    atom_out["x"] = mean (list_x)
    atom_out["y"] = mean (list_y)
    atom_out["z"] = mean (list_z)
    
    return atom_out
    

def countAtom (dico_count, PDB_parsed, debug = 0):
    """
    Retrieve number of acceptor and donor hydrogene bond
    args: -> dico count
          -> Pocket parsed
    return: -> integers count
    """
    count = 0
    
    for atom in PDB_parsed : 
        residue = tool.transformAA(atom["resName"])
        if debug : print residue
        
        if residue in dico_count : 
            atom_Name = atom["name"]
            if atom_Name in dico_count[residue] : 
                count = count + 1
    return count

    

def hydrophobicityKyte (flag, seq, filout, debug = 1):
    """Calcul hydrophobicity index
    args: -> flag
          -> seq pocket
          -> filout
    return: -> NONE write in filout
    """
    nb_res = float(len(seq))
    h_index = 0
    
    for aa in seq:
        h_index = h_index + kyte_hydropathy_index[aa]
    
    if debug : print h_index,nb_res, "Hydrophobic Kyte"
    try : filout.write(flag+"_hydrophobic_kyte\t"+str(h_index/nb_res)+"\n")
    except : filout.write(flag+"_hydrophobic_kyte\tNA\n")


def charge(flag,seq,fileout,na= 0):# a recoder pour etre comme perola qui inclus astartic + glutamic + lysine + arginine + metaux
    if na==0:
        nb_pos = 0
        for aa_pos in res_pos:
            nb_pos += string.count(seq,aa_pos)
        nb_neg = 0
        for aa_neg in res_neg:
            nb_neg += string.count(seq,aa_neg)
        nb_tot = nb_pos-nb_neg
    elif na==1:
        nb_tot = "NA"
    fileout.write(flag+"_charge\t"+str(nb_tot)+"\n")


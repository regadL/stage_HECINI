"""
BORREL Alexandre
04-2012
"""


def listDescriptorFpocket ():
    """
    List descriptors Fpocket
    args: NONE
    return: list of descriptors
    """
    
    # list_out = ["Pocket Score", "Number of V. Vertices", "Mean alpha-sphere radius", "Mean alpha-sphere SA", "Mean B-factor", "Hydrophobicity Score", "Polarity Score", "Volume Score", "Real volume", "Charge Score", "Local hydrophobic density Score", "Number of apolar alpha sphere", "Proportion of apolar alpha sphere"]
#    list_out = ["Number of V. Vertices", "Mean alpha-sphere radius", "Mean alpha-sphere SA", "Mean B-factor", "Real volume", "Number of apolar alpha sphere", "Proportion of apolar alpha sphere"]
    list_out = ["Mean alpha-sphere radius", "Mean alpha-sphere SA", "Mean B-factor", "Real volume", "Proportion of apolar alpha sphere"]
    return list_out


def dictionaryDescriptorFpocket ():
    """Generate stock structure
    arg: NULL
    return: dictionary"""
    
    dico_out = {}
    
# File atm.pdb header
    dico_out ["Pocket Score"] = []
    dico_out ["Drug Score"] = []
    dico_out ["Number of V. Vertices"] = []
    dico_out ["Mean alpha-sphere radius"] = []
    dico_out ["Mean alpha-sphere SA"] = []
    dico_out ["Mean B-factor"] = []
    dico_out ["Hydrophobicity Score"] = []
    dico_out ["Polarity Score"] = []
    dico_out ["Volume Score"] = []
    dico_out ["Real volume"] = []
    dico_out ["Charge Score"] = []
    dico_out ["Local hydrophobic density Score"] = []
    dico_out ["Number of apolar alpha sphere"] = []
    dico_out ["Proportion of apolar alpha sphere"] = []



# file info Fpocket    
#    dico_out ["Score"] = []
#    dico_out ["Druggability Score"] = []
#    dico_out ["Number of Alpha Spheres"] = []
#    dico_out ["Total SASA"] = []
#    dico_out ["Polar SASA"] = []
#    dico_out ["Apolar SASA"] = []
#    dico_out ["Volume"] = []
#    dico_out ["Mean local hydrophobic density"] = []
#    dico_out ["Mean alpha sphere radius"] = []
#    dico_out ["Mean alp. sph. solvent access"] = []
#    dico_out ["Apolar alpha sphere proportion"] = []
#    dico_out ["Hydrophobicity score"] = []
#    dico_out ["Volume score"] = []
#    dico_out ["Polarity score"] = []
#    dico_out ["Charge score"] = []
#    dico_out ["Proportion of polar atoms"] = []
#    dico_out ["Alpha sphere density"] = []
#    dico_out ["Cent. of mass - Alpha Sphere max dist"] = []
#    dico_out ["Flexibility"] = []
    
    return dico_out

def dictionnaryDescriptorCalcul ():
    """NOT USE -> DELL"""
    
    dico_out= {}
    # Area descriptors
    dico_out ["planarity"] = []
    dico_out ["d1+d2"] = []
    dico_out ["d4+d5"] = []
    dico_out ["narrowness"] = []
    dico_out ["rugosity"] = []
    # Atomic descriptors
    dico_out ["residues"] = []
    dico_out ["A"] = []
    dico_out ["C"] = []
    dico_out ["D"] = []
    dico_out ["E"] = []
    dico_out ["F"] = []   
    dico_out ["G"] = []
    dico_out ["H"] = [] 
    dico_out ["I"] = [] 
    dico_out ["K"] = []   
    dico_out ["L"] = []
    dico_out ["M"] = [] 
    dico_out ["N"] = [] 
    dico_out ["P"] = [] 
    dico_out ["Q"] = []   
    dico_out ["R"] = []
    dico_out ["S"] = [] 
    dico_out ["T"] = [] 
    dico_out ["V"] = [] 
    dico_out ["W"] = [] 
    dico_out ["Y"] = []   
    dico_out ["aromatic_residues"] = [] 
    dico_out ["polar_residues"] = []  
    dico_out ["tiny_residues"] = []
    dico_out ["hydrophobic_residues"] = []    
    dico_out ["aliphatic_residues"] = []   
    dico_out ["negative_residues"] = [] 
    dico_out ["positive_residues"] = [] 
    dico_out ["charged_residues"] = []
    dico_out ["carbone"] = []
    dico_out ["carbone_mainchain"] = [] 
    dico_out ["carbone_sidechain"] = []  
    dico_out ["azote"] = []
    dico_out ["azote_mainchain"] = [] 
    dico_out ["azote_sidechain"] = []  
    dico_out ["soufre"] = []
    dico_out ["soufre_mainchain"] = [] 
    dico_out ["soufre_sidechain"] = []  
    dico_out ["oxygene"] = []
    dico_out ["oxygene_mainchain"] = [] 
    dico_out ["oxygene_sidechain"] = [] 
    dico_out ["hydrogene"] = []
    dico_out ["hydrogene_mainchain"] = [] 
    dico_out ["hydrogene_sidechain"] = [] 
    # energy descriptors
    dico_out ["hbond_acceptor"] = [] 
    dico_out ["hbond_donor"] = []
    dico_out ["charge"] = [] 
    dico_out ["polarity_ratio"] = [] 
    dico_out ["hydrophobicity_ratio"] = [] 
    # volume descriptors
    dico_out ["volume"] = [] 
    dico_out ["surface"] = []
    dico_out ["compacity"] = []     
    dico_out ["lambda0"] = [] 
    dico_out ["lambda1"] = []
    dico_out ["lambda2"] = []
    dico_out ["longueur0"] = [] 
    dico_out ["longueur1"] = []
    dico_out ["longueur2"] = []    
    
    return dico_out
    


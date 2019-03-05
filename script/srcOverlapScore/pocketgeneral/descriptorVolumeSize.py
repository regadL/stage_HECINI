"""
BORREL Alexandre
04-2012
"""

from os import system, path
from re import search
from math import sqrt, pi
from numpy import mean, var, array, dot, linalg, transpose, arctan


import runOtherProg
from PDB import *


def volumeChimera(filout, path_file_protomol_mol2, path_dir_result):
    """
    Calcul chimera descriptor with protomol (mol2)
    args: -> filout descriptors
          -> path file protomol (mol2)
          -> path directory result
    return: NONE write file descriptor (filout)
    """
    
    # write script for chimera
    path_script_chimera = writeChimeraScript (path_file_protomol_mol2, path_dir_result)
    path_result_chimera = runOtherProg.chimeraVolume (path_script_chimera, path_dir_result)
    
    result_parsed = parseChimeraResult (path_result_chimera)
    
    for descriptor in result_parsed.keys () :
        filout.write ("pocket_" + descriptor + "\t" + str (result_parsed[descriptor]) + "\n")
    
#    system ("rm " + path_result_chimera)
#    system ("rm " + path_script_chimera)
    
def writeChimeraScript (path_file_protomol, path_dir_result) : 
    """
    Write script for Chimera
    args: -> path file protomol
          -> path directory result
    return: -> path script chimera
    """
    
    path_script_chimera = path_dir_result + "script_chimera"
    file_script = open (path_script_chimera, "w")
    file_script.write ("\nopen " + path_file_protomol + "\nsurf #0\nmeasure volume #0\nmeasure area #0\n\n\n\n")
    file_script.close()
    return path_script_chimera


def parseChimeraResult (path_result_chimera) : 
    
    # initialization out put
    dico_out = {}
    dico_out["volume"] = "NA"
    dico_out["surface"] = "NA"
    dico_out["compacity"] = "NA"
    
    if path.getsize(path_result_chimera) == 0 : 
        return dico_out
    
    filin = open (path_result_chimera, "r")
    list_line_chimera = filin.readlines ()
    filin.close ()
    
    for line_chimera in list_line_chimera : 
        if search("^MSMS", line_chimera) : 
            if search ("volume", line_chimera) : 
                dico_out["volume"] = float(line_chimera.split ("=")[-1].strip ().replace (" ", ""))
            elif search ("area", line_chimera) : 
                dico_out["surface"] = float(line_chimera.split ("=")[-1].strip ().replace (" ", ""))
    
    try :dico_out["compacity"] =  dico_out["volume"] / dico_out["surface"]
    except : pass
    
    return dico_out




def compute_lambdas(path_file_pocket, path_dir_pocket, path_file_protomol_pdb, path_pdb, filout, debug = 0):
        """
        Surement a reprendre car le commentaire met sur le ligand et on est sur la pochet ????
        caculer sur le protomol ne comprend ce que fait la premier partie du script sur la poche
        car non utilise ensuite je pense bout de code perdu dans les programmes ou subtitue entre
        ligand et protomol
        """
        # calcul du nombre de sondes dans le fichier protomol ou ligand
        # ATTENTION : il y a probablement des protomol ou des ligands en ATOM !!
        system("awk '$1==\"ATOM\" || $1==\"HETATM\" {print $0}' "+path_file_pocket+" | wc -l > " + path_dir_pocket + "temp_sonde.out")
        fileS = open(path_dir_pocket + "temp_sonde.out",'r')
        nS = fileS.readline()
        nS = string.replace(nS,'>\n','')
        fileS.close()
        system("rm " + path_dir_pocket + "temp_sonde.out")
        nS = int(nS)
        # ouverture du fichier par PDBpy
        protomol = PDB(path_file_protomol_pdb)  ###filein  = ligand.pdb
        x = []
        y = []
        z = []
        # recuperation des coordonnees
        for res in protomol:
            for at in res:
                xyz = at.xyz()
                x.append(xyz[0])
                y.append(xyz[1])
                z.append(xyz[2])
        # equivaut au scale(center=T,scale=F) de R ie matrice centree mais non reduite
        xm = x - mean(x)
        ym = y - mean(y)
        zm = z - mean(z)
        xv = xm/sqrt(var(x))
        yv = ym/sqrt(var(y))
        zv = zm/sqrt(var(z))
        # matrice centree non reduite
        scale = array([xm,ym,zm])
        # matrice centree reduite
        #scale = array([xv,yv,zv])
        # matrice des variances covariances
        covar = dot(scale,transpose(scale))
        lambdas = linalg.eig(covar)[0]
        sum_lambda = lambdas[0] + lambdas[1] + lambdas[2]
        lambda0 = int(lambdas[0])/float(sum_lambda)
        lambda1 = int(lambdas[1])/float(sum_lambda)
        lambda2 = int(lambdas[2])/float(sum_lambda)
        # tri
        sort = [lambda0,lambda1,lambda2]
        sort.sort()
        # mini - maxi pour calculer les descripteurs de Nayal
        #lambdamax = maxi3(lambda0,lambda1,lambda2)
        #lambdamin = mini3(lambda0,lambda1,lambda2)
        lambdamax = sort[2]
        lambdamin = sort[0]
        lambdaint = sort[1]
        filout.write("pocket_lambda0\t%.2f" %lambdamin+"\n")
        filout.write("pocket_lambda1\t%.2f" %lambdaint+"\n")
        filout.write("pocket_lambda2\t%.2f" %lambdamax+"\n")
        
        gcur = lambdamax*lambdamin
        gcur = "%.2f" %gcur
        mcur = (lambdamax+lambdamin)/2.
        mcur = "%.2f" %mcur
        curved = sqrt((lambdamax**2+lambdamin**2)/2.)
        curved = "%.2f" %curved
        index = -2*arctan((lambdamax+lambdamin)/float(lambdamax-lambdamin))/pi
        index = "%.2f" %index
        # de cote pour le moment ... attention dans la lecture finale ! ne pas prendre que le -1
        vps = linalg.eig(covar)[1]
#        vp0 = vps[0]
#        vp1 = vps[1]
#        vp2 = vps[2]
        #filout.write(flag+"_vp0\t"+str("%.2f"%vp0[0])+" "+str("%.2f"%vp0[1])+" "+str("%.2f"%vp0[2])+"\n")
        #filout.write(flag+"_vp1\t"+str("%.2f"%vp1[0])+" "+str("%.2f"%vp1[1])+" "+str("%.2f"%vp1[2])+"\n")
        #filout.write(flag+"_vp2\t"+str("%.2f"%vp2[0])+" "+str("%.2f"%vp2[1])+" "+str("%.2f"%vp2[2])+"\n")
        
        distance = []
        min_var = 10000000.
        for v in range(len(vps)):
            # calcul de la "profondeur" sur vp0
            # le vecteur orthogonal definissant le plan
            n1 = vps[v][0]
            n2 = vps[v][1]
            n3 = vps[v][2]
            # un point quelconque qui sert a definir le plan
            xA =  x[0]
            yA =  y[0]
            zA =  z[0]
            # equation du plan ax+by+cz+d=0 : ensemble des points M tels que n.AM = 0
            # donc en fait n1x+n2y+n3z-(n1xA+n2yA+n3zA) = 0
            # calcul de la distance de toutes les sondes au plan pour trouver la max
            maxD = 0.
            maxI = 0
            for i in range(len(x)):
                D = abs(n1*x[i]+n2*y[i]+n3*z[i]-n1*xA-n2*yA-n3*zA)/sqrt(n1**2+n2**2+n3**2)
                if D > maxD:
                    maxD = D
                    maxI = i
        # on deplace le plan au point maxI (donc a une extremite) et on recherche le max qui est le vrai cette fois
        # et qui correspond a l'elongation max
            xI = x[maxI]
            yI = y[maxI]
            zI = z[maxI]
            profD = 0.
            profI = 0
            for i in range(len(x)):
                D = abs(n1*x[i]+n2*y[i]+n3*z[i]-n1*xI-n2*yI-n3*zI)/sqrt(n1**2+n2**2+n3**2)
                if D > profD:
                    profD = D
                    profI = i
                    #print "distance max au plan :","%.2f"%profD
            distance.append(profD)
        
            # ATTENTION -> retir profondeur
        
        if debug : print distance
        filout.write("pocket_longueur0\t%.2f" %distance[0]+"\n")
        filout.write("pocket_longueur1\t%.2f" %distance[1]+"\n")
        filout.write("pocket_longueur2\t%.2f" %distance[2]+"\n")





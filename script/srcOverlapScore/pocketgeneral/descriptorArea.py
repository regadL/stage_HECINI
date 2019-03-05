"""
BORREL Alexandre
04-2012
calcul descriptor
"""

from math import sqrt, log, pi
from Scientific.Geometry.Objects3D import Vector, Plane, Line
from Scientific.Functions.LeastSquares import leastSquaresFit
from numpy import median

from PDB import *
import runOtherProg
import parseNACCESS
import parsePDB


def narrowness(filein,fileout,na=0):
    """
    Calcul narrowness
    """
    if na==0:
        params = least_square_plane(filein)
        A = params[0]
        B = params[1]
        C = params[2]
        #div = A**2+B**2+1
        p1 = Vector(0.,0.,C)
        p2 = Vector(0.,1.,B+C)
        p3 = Vector(1.,0.,A+C)
        plan = Plane(p1,p2,p3)
        pdb = PDB(filein)
        proj = []
        #proj1 = []
        for res in pdb:
            for at in res:
                xyz = at.xyz()
                x = xyz[0]
                y = xyz[1]
                z = xyz[2]
                point = Vector(x,y,z)
                coordproj = Plane.projectionOf(plan,point)
                #xproj = (x*(1+B**2)-A*B*y+A*z-A*C)/div
                #yproj = (-A*B*x+(A**2+1)*y+B*z-B*C)/div
                #zproj = (A*x+B*y+(A**2+B**2)*z+C)/div
                proj.append(coordproj)
        dmax = 0
        imax = 0
        jmax = 0
        for i in range(len(proj)-1):
            for j in range(i+1,len(proj)):
                d = sqrt((proj[i][0]-proj[j][0])**2+(proj[i][1]-proj[j][1])**2+(proj[i][2]-proj[j][2])**2)
                if d > dmax:
                    dmax = d
                    imax = i
                    jmax = j
        d3 = dmax
        pointI = Vector(proj[imax][0],proj[imax][1],proj[imax][2])
        pointJ = Vector(proj[jmax][0],proj[jmax][1],proj[jmax][2])
        direction = pointJ-pointI
        droite = Line(pointI,direction)
        #vecteurs_ortho = []
        vecteurs_dist = []
        dist_droite = []
        dist_gauche = []
        pointM = proj[0]
        projM = droite.projectionOf(pointM)
        for p in proj:
            x = p[0]
            y = p[1]
            z = p[2]
            point = Vector(x,y,z)
            projete = droite.projectionOf(point)
            D = droite.distanceFrom(point)
            scalaire = (point-projete)*(pointM-projM)
            vecteurs_dist.append(D)
            if scalaire > 0:
                dist_gauche.append(D)
            else:
                dist_droite.append(D)
        # a revoir : pb si d5 vide ... exemple : 1nt4 -> write NA if list -> empty
        if dist_droite == [] or dist_gauche == [] : 
            fileout.write ("pocket_narrowness\tNA\n")
            fileout.write("pocket_d4+d5\tNA\n")
            return 
        
        d4 = max(dist_droite)
        d5 = max(dist_gauche)
        d4d5 = d4+d5
        d4d5 = "%.2f" %d4d5
        narrow = 1-((d4+d5)/d3)
        narrow = "%.2f" %narrow
        #print narrow
    elif na==1:
        narrow = "NA"
        d4d5 = "NA"
    fileout.write("pocket_narrowness\t"+str(narrow)+"\n")
    fileout.write("pocket_d4+d5\t"+str(d4d5)+"\n") # for narrowness calculation


def planarity(file_protomol,fileout,na=0):
    """
    calcul planarity fonction file protomol
    calcul difference betwen plane and other cote
    
    """
    if na==0:
        params = least_square_plane(file_protomol)
        A = params[0]
        B = params[1]
        C = params[2]
        pdb = PDB(file_protomol)
        dist = []
        distij = []
        for res1 in pdb:
            for at1 in res1:
                xyz1 = at1.xyz()
                x1 = xyz1[0]
                y1 = xyz1[1]
                z1 = xyz1[2]
                d = (A*x1+B*y1-z1+C)/sqrt(A**2+B**2+1)
                dist.append(d)
                for res2 in pdb:
                    for at2 in res2:
                        xyz2 = at2.xyz()
                        x2 = xyz2[0]
                        y2 = xyz2[1]
                        z2 = xyz2[2]
                        dij = sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
                        distij.append(dij)
        d1 = abs(min(dist))
        d2 = max(dist)
        d1d2 = d1+d2
        d1d2 = "%.2f" %d1d2
        plan = 1-((d1+d2)/max(distij))
        plan = "%.2f" %plan
        #print plan
    elif na==1:
        plan = "NA"
        d1d2 = "NA"
    fileout.write("pocket_planarity\t"+str(plan)+"\n")
    #fileout.write("pocket_d1+d2\t"+str(d1d2)+"\n") # directly correlated with size of pocket

def least_square_plane(filein):
    pdb = PDB(filein)
    data = []
    for res in pdb:
        for at in res:
            xyz = at.xyz()
            x = xyz[0]
            y = xyz[1]
            z = xyz[2]
            data.append(((x,y),z))
    #print data
    parameters = [data[1][0][0],data[1][0][1],data[1][1]]
    params, chisquared = leastSquaresFit(model,parameters,data)
    return params


def model(params,xy):
    c0, c1, c2 = params
    x, y = xy
    return c0*x + c1*y + c2


def rugosity(path_file_pocket_atom, path_complex, path_dir_result, PDB_ID, filout, debug = 0):
    """
    Calcul rugosity Pettit et al 1999, based neighbor atoms that exposed solvant atoms 
    THREHOLD -> 5A fix in publication
    args: -> path file pocket residue acc
          -> filout for write descriptors
          -> path directory dataset
          -> path directory result
          -> PDB ID
    return: NONE write in filout
    """
    list_safd = []
    list_probes = [1.0,1.1]
    pathes_files_asa = {}
    for probe in list_probes : 
        pathes_files_asa[probe] = runOtherProg.runNACESSWithOption(path_complex, probe, path_dir_result)
    if debug : print pathes_files_asa
    
    complex_parsed = parsePDB.loadCoordSectionPDB(path_complex)
    pocket_parsed = parsePDB.loadCoordSectionPDB(path_file_pocket_atom)
    appendAccParsed (complex_parsed, pathes_files_asa)
    appendAccParsed (pocket_parsed, pathes_files_asa)
    
    
    for atom_pocket in pocket_parsed : 
        if debug : 
            for dico_keys in atom_pocket.keys () : 
                print dico_keys, atom_pocket[dico_keys]
        contact_area = []
        for probe in list_probes : 
            sum_ABS = 0.0
            if atom_pocket[probe]["ABS"] > 0.0 :
                for atom_complex in complex_parsed :
                    if parsePDB.distanceTwoatoms(atom_complex, atom_pocket) < 5 :
                        sum_ABS = sum_ABS + atom_complex[probe]["ABS"]
            contact_area.append("%.2f" %sum_ABS)
        try : 
            safd = calculSAFD(contact_area,list_probes)
            list_safd.append(safd)
        except :
            pass
        
    median_safd = "%.2f" %median(list_safd)
    if debug : print median_safd, "Rugosity"    
    filout.write("pocket_rugosity\t" + median_safd + "\n")   
     



def appendAccParsed (complex_parsed, pathes_files_asa) : 
    """
    Append in dictionary PDB parsed value of ASA by atoms
    args: -> complex parsed with parsePDB parser
          -> dictionary with pathes file asa
    return: NONE change complex parsed
    """
    for probe in pathes_files_asa.keys () :
        asa_parsed = parseNACCESS.fileASA(pathes_files_asa[probe])
        for atom_complex in complex_parsed :
            nb_atom_asa = len (asa_parsed)
            i = 0 
            while i < nb_atom_asa : 
                if atom_complex["serial"] == asa_parsed[i]["atomSeq"] and atom_complex["chainID"] == asa_parsed[i]["chainID"] and atom_complex["name"] == asa_parsed[i]["name"]:
                    atom_complex[probe] = {}
                    atom_complex[probe]["ABS"] = asa_parsed[i]["ABS"]
                    del asa_parsed[i]
                    break
                else :
                    i = i + 1
                    
        # check if probe append in parsing
        # Acid nuc not present in asa file and run pqr bug
        nb_atom = len (complex_parsed)
        i = 0
        while i< nb_atom : 
            if not probe in complex_parsed[i].keys () : 
                del complex_parsed[i]
                nb_atom = nb_atom - 1
            else :
                i = i + 1
            
        
    
def calculSAFD(vecteur,sonde):
    a = (log(float(vecteur[1]))-log(float(vecteur[0])))/(log(float(sonde[1]))-log(float(sonde[0])))
    b = log(float(vecteur[0]))-a*log(float(sonde[0]))
    return(2-a)


def atomsAround(atomei,path_file_pdb, extend = 5):
    """
    retrieve atom closing
    """
    
    liste_atomes = []
    prot = PDB(path_file_pdb)
    for res in prot:
        for at in res:
            if int(at.atmNum()) == int(atomei):
                xyz = at.xyz()
                xi = xyz[0]
                yi = xyz[1]
                zi = xyz[2]
                break
    for res in prot:
        for at in res:
            xyzc = at.xyz()
            xc = xyzc[0]
            yc = xyzc[1]
            zc = xyzc[2]
            dist = sqrt((xc-xi)**2+(yc-yi)**2+(zc-zi)**2)
            a = int(at.atmNum())
            if at.header() == 'ATOM':
                if (a != int(atomei) and dist <= extend):
                    liste_atomes.append(a)
    return liste_atomes
---
title: "Etude des poches de PR2"
author: "Leslie REGAD et Akram HECINI"
date: '2019-03-06'
output:
  html_document:
    code_folding: show #hide
    self_contained: yes
    fig_caption: yes
    highlight: pygments #pour les sorties R
    theme: spacelab
    toc: yes  #for add table of contents
    toc_depth: 3
    toc_float: yes
    df_print: paged
    keep_md: yes  #pour garder le .md after run
  slidy_presentation:
    smart: no
    slide_level: 2
    self_contained: yes
    fig_caption: yes
    fig_height: 6
    fig_width: 7
    highlight: tango
    incremental: no
    keep_md: yes
    smaller: yes
    theme: cerulean
    toc: yes
    widescreen: yes
  ioslides_presentation:
    slide_level: 2
    self_contained: no
    colortheme: dolphin
    fig_caption: yes
    fig_height: 6
    fig_width: 7
    fonttheme: structurebold
    highlight: tango
    smaller: yes
    toc: yes
    widescreen: yes
  pdf_document:
    fig_caption: yes
    highlight: zenburn
    toc: yes
    toc_depth: 3
  beamer_presentation:
    colortheme: dolphin
    fig_caption: yes
    fig_height: 6
    fig_width: 7
    fonttheme: structurebold
    highlight: tango
    incremental: no
    keep_tex: no
    slide_level: 2
    theme: Montpellier
    toc: yes
font-import: http://fonts.googleapis.com/css?family=Risque
font-family: Garamond
transition: linear
---




# Objectif du projet

Le but du projet est d'étudier les différents sites de liaison de la protéase du VIH-2 (PR2) dans le but d'identifier de nouveau site de liaison



# Data disponibles 

* le répertoire `data/PR2_19` contient les fichiers PDB des 19 structures de PR2 disponibles dans la PDB.


# Protocole

* Estimation des poches des 19 PR2 à l'aide de Fpocket  
    + étape 1 : A l'aide du programme Fpocket, identifier toutes les poches des 19 PR2. 
    + étape 2 : calcul du score de druggabilité de chaque poche à l'aide du logiciel PockDrug.  
    + étape 3 : calcul du chevauchement entre toutes les poches, pour identifier les poches similaires des différentes protéines.  

* Estimation des poches des 19 PR2 à l'aide de FTMap  
    + étape 1 : A l'aide du programme FTMap, identifier toutes les poches des 19 PR2. 
    + étape 2 : comparer les résultats avec les poches obtenues avec Fpocket
    
* Classification des poches basées sur leur description  
    + étape 1 : calculer les descripteurs de différentes poches estimées
    + étape 2 : classification des poches

* Selection des poches d'intérêt

* Etude des interactions entre les résidus de ces poches et les autres résidus des PR2.


#04/03/2019 


# Estimation des poches

- faite en utilisant le programme fpocket.

- Génération du programme `/home/hecini/Research/stage_HECINI/data/PR2_19/fpocket.sh` pour lancer l'estimation des poches sur les 19 PR2.

- les fichiers output se trouvent dans :  `/home/hecini/Research/stage_HECINI/data/PR2_19/pockets`

- Génération du programme copier_coller.sh `/home/hecini/Research/stage_HECINI/data/PR2_19/pockets/copier_coller.sh`  qui permet mettre les pockets et le fichier pdb en question dans un seul dossier. 

- Visualisation des pockets via pymol. pour chaque protèine il y'a une session *.pse qui a été sauvegardé afin de revoir les résultats facilement. 

# Dénombrement des poches par protéine

fait manuellement

[résultats](pdf/nombre_poches.pdf)

pdb | poche | poche prin |
----|-------|------------|
1HSI|   2   |    0       |
1HSI|   2   |    0       |
1HII| 5     |            |
1JLD|10     |            |
1IDB|9      |            | 
2HPE|11 ?   |            |
3EC0|9 0    |            |
4UPJ|6      | 0 et 1     |
1HSH|7      |0 et 2     |
1IVP|9 
3ECG| 11| 0|
5UPJ| 6| 0|
1HSI| 10| //|
1IVQ |7| 0|
2MIP |7 ||
3S45 |10 |0 et 1
6UPJ |6 | 0 1 et 8
1IDA |5 | 0 1 et 2
3EBZ |10 |0|
3UPJ |7 |0

# Le 05/03/2019 et le 06/03/19 

#Création d'un nouveau répertoir contenant deux sous_repertoires 

    --> Pockets qui contient les fichiers pdb des poches. Les fichier pdb ont été renomés en utilisant 
        un script bash rename.sh
        
    --> proteins : qui contient les 19 strcutures renomées avec rename.sh
    
    -->rename.sh est modifié selon le fichier que l'on souhaite renomer. 
    
  
#Calcul des scores et analyse de données 

- Développer un programme python pour calculer les scores de similarité entre les poches
    [score_v2.py](script/score_v2.py) 
- génération d'un fichier [scores.csv](script/score.csv) 
- Préparation du fichier csv  avec R pour analyser les données:

    --> diviser la 1ere colonne en deux pour séparer les noms des poches
    
    --> Transformation du dataFrame en Matrix numérique 
    
    --> générer un pheatmap
  
  En utilisant [ ce script.R ](script/scores_analysis.R) 
  

![pheatmap des poches](script/Rplot.png)

  
  
#Bibliotheque :

  lecture de "Damm KL, Ung PM, Quintero JJ, Gestwicki JE, Carlson HA. A poke in the eye: inhibiting HIV-1 protease through its flap-recognition pocket. Biopolymers. 200889:643-52"



  
  


```r
pdb <- c("1HSI", "3S45")
nbr.pocket <- c(10,15)
```

- etudier la repartition du nombre de poches par protéine.



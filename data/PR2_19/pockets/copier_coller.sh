#!/bin/bash
#Mettre le fichier out.pdb avec l'ensemble des pockets pour la visualisation avec pymol 

for file in *_out
do
	cd $file
	cp *.pdb pockets
	cd .. 

done

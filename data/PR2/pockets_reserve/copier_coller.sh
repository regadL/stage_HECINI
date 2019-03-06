#!/bin/bash
#Mettre le fichier out.pdb avec l'ensemble des pockets pour la visualisation avec pymol 

for file in *_out/
do
	cd $file
		for subfile in $file
			cd subfile 

			cp *.pdb /home/hecini/Research/stage_HECINI/data/PR2/pockets
			cd ..
	cd ..

done

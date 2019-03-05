#!/usr/bin/env Rscript


# By BORREL Alexandre
# 04-2012
require(graphics)
source ("tool.R")


generateColor= function(nb_color){

	if (nb_color == 2){
		return (c("#66FF77","#FF0000"))
	}
	else if(nb_color == 4){
		return (c("#66FF77","#338800","#FFCCCC", "#FF0000"))
	}
	else if (nb_color == 6){
		return (c("#003300","#338800","#66FF77","#FFFFFF","#FFCCCC", "#FF0000","#8B2323"))
	}
	else if (nb_color == 10){
		return (c("#6600CC","#0000CC","#3366FF","#33CCFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFFFF","#FFFF66","#FFCC33","#FF6633","#CC0000"))

	}
}


generateLegend = function (length_tot, cut){
	
	list_out = seq(0,length_tot, cut)
	if (list_out[length(list_out)] != length_tot){
		list_out = append(list_out, length_tot)
	}
	return (list_out)
}


generatePosition = function (list_legend, value_ecart){

	list_out = c()
	for (element in list_legend){
		list_out = append(list_out,element * value_ecart)
	}
	return (list_out)
}




cardMatrixCor = function(matrixIN, name_file, nb_color){

	nb_col = dim(matrixIN)[2]
	nb_line = dim(matrixIN)[1]

	y_names = seq(0,nb_col,1)
	x_names =  seq(0,nb_line,1)
	

	list_color = generateColor(nb_color)
	
	#postscript (file = paste (name_file, ".ps", sep = ""), width = 0, height = 0)

	#par( mar=c(25,25,1.5,1.5))
	#image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color)
	#grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
	#box()
	# place les petites barres 
	#axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
	#axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)

	# place les positions en fonction du cut
	# ecart1 = 1/(nb_line-1)	
	# ecart2 = 1/(nb_col-1)
	# list_L1 = generateLegend (nb_line,1)
	# list_L2 = generateLegend (nb_col,1)
	
	# place les legendes
	# posX = generatePosition(list_L1, ecart1)
	# posY = generatePosition(list_L2, ecart2)
	# axis(1,seq(0,1,(1/(nb_line-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
	# axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
	# dev.off()




	# PNG file

	# list descriptor

	# fix breaks
	x = 0
	while ((length (list_color)+1) != length (seq(-1,1, 2/(nb_color+x)))){
		x = x + 1
	}


	png (file = paste (name_file, ".png", sep = ""), width = 1600, height = 1600)

	par( mar=c(25,25,1.5,1.5))
	image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color, breaks = seq(-1,1, 2/(nb_color+x)))
	grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
	box()
	# place les petites barres 
	axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
	axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)

	# place les positions en fonction du cut
	ecart1 = 1/(nb_line-1)	
	ecart2 = 1/(nb_col-1)
	list_L1 = generateLegend (nb_line,1)
	list_L2 = generateLegend (nb_col,1)
	
	# place les legendes
	posX = generatePosition(list_L1, ecart1)
	posY = generatePosition(list_L2, ecart2)
	axis(1,seq(0,1,(1/(nb_line-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
	axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
	dev.off()

}



cardMatrix = function(matrixIN, name_file, nb_color){

	nb_col = dim(matrixIN)[2]
	nb_line = dim(matrixIN)[1]

	list_color = generateColor(nb_color)
	
	#postscript (file = paste (name_file, ".ps", sep = ""), width = 0, height = 0)

	#par( mar=c(25,25,1.5,1.5))
	#image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color)
	#grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
	#box()
	# place les petites barres 
	#axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
	#axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)

	# place les positions en fonction du cut
	# ecart1 = 1/(nb_line-1)	
	# ecart2 = 1/(nb_col-1)
	# list_L1 = generateLegend (nb_line,1)
	# list_L2 = generateLegend (nb_col,1)
	
	# place les legendes
	# posX = generatePosition(list_L1, ecart1)
	# posY = generatePosition(list_L2, ecart2)
	# axis(1,seq(0,1,(1/(nb_line-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
	# axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
	# dev.off()




	# PNG file

	# list descriptor
	#list_desc = c("hydrophobic_kyte", "charge", "hydrophobicity_ratio_pocket", "polarity_ratio_pocket", "PSI", "RADIUS_CYLINDER", "X._ATOM_CONVEXE", "CONVEX.SHAPE_COEFFICIENT", "INERTIA_2", "INERTIA_3", "INERTIA_1", "HEIGHT_MIN_CYLINDER", "FACE", "HEIGHT_CYLINDER", "PCI", "DIAMETER_HULL", "SURFACE_HULL", "RADIUS_MIN_CYLINDER", "VOLUME_HULL", "SMALLEST_SIZE", "RADIUS_HULL" , "c_residues", "c_atom", "p_charged_residues", "p_positive_residues", "p_pro_residues", "p_polar_residues", "p_tiny_residues", "p_aliphatic_residues", "p_negative_residues", "p_aromatic_residues", "p_hydrophobic_residues", "p_ND1_atom", "p_hbond_acceptor_atom", "p_Nlys_atom","p_nitrogen_atom", "p_Cgln_atom", "p_Car_atom", "p_N_atom", "p_hyd_atom", "p_hbond_donor_atom", "p_sulfur_atom", "p_Ccoo_atom", "p_NE2_atom", "p_side_chain_atom",  "p_Otyr_atom", "p_O_atom", "p_Carg_atom", "p_Ooh_atom", "p_Ntrp_atom", "p_hydrophobic_atom", "p_main_chain_atom", "p_oxygen_atom", "P_C_atom", "p_carbone_atom", "p_Ocoo_atom", "p_S_atom") 
	#print (length (list_desc))

	png (file = paste (name_file, ".png", sep = ""), width = 1600, height = 1600)

	par( mar=c(25,25,1.5,1.5))
	image(as.matrix(matrixIN), yaxt = "n", xaxt = "n", col = list_color)
	grid(nx = nb_line, ny = nb_col, col = "black", lwd = 1, lty = 1)
	box()
	# place les petites barres 
	axis(1,seq(0,1,(1/(nb_line-1))), labels = FALSE)
	axis(2,seq(0,1,(1/(nb_col-1))), labels = FALSE)

	# place les positions en fonction du cut
	ecart1 = 1/(nb_line-1)	
	ecart2 = 1/(nb_col-1)
	list_L1 = generateLegend (nb_line,1)
	list_L2 = generateLegend (nb_col,1)
	
	# place les legendes
	posX = generatePosition(list_L1, ecart1)
	posY = generatePosition(list_L2, ecart2)
	axis(1,seq(0,1,(1/(nb_line-1))),rownames (matrixIN), cex.axis = 1.75, las = 2)
	axis(2,seq(0,1,(1/(nb_col-1))),colnames (matrixIN), cex.axis = 1.75, las = 2)
	dev.off()

}





is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}

#########
# MAIN  #
#########

args = commandArgs(TRUE)

path_matrix = args[1]
nb_color = as.integer(args[2])
cor_option = as.integer(args[3])

data = read.table (path_matrix, header = T, sep = "\t")

# gestion du CSV
#data = read.csv (path_matrix,header = T, sep = "," )
#rownames (data) = data[,1]
#data = data[,-1]


data = as.matrix(data)

data = try(na.omit(data)) # deleted line with NA

if (is.integer0(which(colnames(data) == "drugg")) == FALSE){
	data = data[,-which (colnames(data) == "drugg")]
}

# order descriptor
d.col = colorACPByTypeOfDescriptors()
data = t(orderByType (d.col, t(data)))
print (dim (data)[1] == dim (data)[2])

if (cor_option == 1){

	# control SD = 0
	sd_desc = apply (data[,1:(dim(data)[2])], 2, sd)
	if (is.integer0(which(sd_desc == 0)) == FALSE){
		data = data[,-which(sd_desc == 0)]
	}	
	data_cor = cor (data)
	cardMatrixCor (data_cor, paste(path_matrix, "_cor", sep = ""), nb_color)
}else {
	
	cardMatrix (data, path_matrix, nb_color)
}

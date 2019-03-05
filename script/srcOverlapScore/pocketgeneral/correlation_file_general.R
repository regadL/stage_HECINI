#!/usr/bin/env Rscript


suppDiff = function (data, list_union){
	data_temp = NULL
	for (PDB in list_union){
		data_temp = rbind (data_temp, data[which (rownames(data) == PDB),])
	}

	return (data_temp) 
}



#########
# MAIN  #
#########

args <- commandArgs(TRUE)
path_file1 = args[1]
path_file2 = args[2]
prefix_name = args[3]



data1 = read.table (path_file1, sep = "\t", header = TRUE)
data2 = read.table (path_file2, sep = "\t", header = TRUE)

data1 = data1[,-which (colnames(data1) == "rugosity")]
data2 = data2[,-which (colnames(data2) == "rugosity")]

list_PDB_1 = rownames (data1)
list_PDB_2 = rownames (data2)
intersect_PDB = intersect(list_PDB_1, list_PDB_2)

data1 = suppDiff (data1, intersect_PDB)
data2 = suppDiff (data2, intersect_PDB)

vect_color = data1$drugg + 1
data1 = data1[,-which (colnames(data1) == "drugg")]
data2 = data2[,-which (colnames(data2) == "drugg")]



list_descriptor = intersect(colnames (data1), colnames (data2))

temp = 1
graph = 1

png (paste(prefix_name, "_", graph, ".png", sep = ""), 2480, 3508)
par (mfrow = c (5, 3))

for (descriptor in list_descriptor){
	if (temp == 16){
		dev.off ()
		graph = graph + 1
		png (paste(prefix_name, "_", graph, ".png", sep = ""), 2480, 3508)
		par (mfrow = c (5, 3))
		temp = 1
	}
	cor_linear = cor.test (data1[,which (colnames(data1) == descriptor)], data2[,which (colnames(data2) == descriptor)])
	plot (data1[,which (colnames(data1) == descriptor)], data2[,which (colnames(data2) == descriptor)], main = as.character(cor_linear$estimate), xlab = descriptor, ylab = descriptor, cex.lab = 2.5, cex.main = 2.5, type = "n")
	text (data1[,which (colnames(data1) == descriptor)], data2[,which (colnames(data2) == descriptor)], labels = rownames (data1), col = vect_color )
	temp = temp + 1 
}
dev.off()

for (descriptor in list_descriptor){
	png (paste(prefix_name, "_", descriptor, ".png", sep = ""), 2400, 1400)
	#par(mar=c(25,25,25,25), oma=c(0,0,0,0))
	cor_linear = cor.test (data1[,which (colnames(data1) == descriptor)], data2[,which (colnames(data2) == descriptor)])
	if (descriptor == "volume"){
		plot (data1[,which (colnames(data1) == descriptor)], data2[,which (colnames(data2) == descriptor)], xlab = "", ylab = "", cex.axis = 3, pch = 15, cex = 2, col = vect_color)
	}
	else if (descriptor == "charged_residues"){
		plot (data1[,which (colnames(data1) == descriptor)], data2[,which (colnames(data2) == descriptor)], xlab = "", ylab = "", cex.axis = 3, pch = 15, cex = 2, col = vect_color )
	}
	else {
		plot (data1[,which (colnames(data1) == descriptor)], data2[,which (colnames(data2) == descriptor)], xlab = "Avec ligand", ylab = "Fpocket", cex.axis = 2.5, cex.lab = 3, pch = 15, cex = 2, type = "n")
		text (data1[,which (colnames(data1) == descriptor)], data2[,which (colnames(data2) == descriptor)], labels = rownames (data1), cex = 2, col = vect_color )	
	}	
	dev.off()
}

warnings()

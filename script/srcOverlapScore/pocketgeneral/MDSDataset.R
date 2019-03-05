#!/usr/bin/env Rscript

# By BORREL Alexandre
# 10-2013

source ("tool.R")


MDSCor = function (data1, data2, path_result, color_point, color_point2, matrix_correspondance){
	#print (data)	

	color_point2[which(color_point2 == 1)] = "#01B0F0"
	color_point2[which(color_point2 == 0)] = "#12A003"

	# fusion matrix
	data = rbind (data1,data2)
	color_point = append (color_point, color_point2)

	#Scale
	data = scale (t(data))

	# correlation
	data.cor=cor(data)

	# dist
	dist1 = 1-data.cor
	#print (dim (dist1))

	# fit MDS
	fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)

	png (paste (path_result, "_corPDB_MDS.png", sep = ""), 2480, 1700)
	par( mar=c(12,12,12,12))
	plot (fit$points[,1], fit$points[,2], main="MDS PDB", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, type = "n")
	text (fit$points[,1], fit$points[,2], labels = colnames (dist1), col = as.integer(color_point),  cex = 2.5)
	
	dev.off()

}


MDSDist = function (data1, data2, path_result, color_point, color_point2, matrix_correspondance){
	#print (data)	

	color_point2[which(color_point2 == 1)] = "#01B0F0"
	color_point2[which(color_point2 == 0)] = "#12A003"

	# fusion matrix
	data = rbind (data1,data2)
	color_point = append (color_point, color_point2)

	#Scale
	ds = scale (data)

	# distance
	dist1=dist(ds)
	

	# fit MDS
	fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
	
	###########################
	# color by type estimator #
	###########################
	color_ALSL = color_point

	color_ALSL[which(color_ALSL == 1)] = 3
	color_ALSL[which(color_ALSL == 2)] = 3
	color_ALSL[which(color_ALSL == "#01B0F0")] = 4
	color_ALSL[which(color_ALSL == "#12A003")] = 4
	

	png (paste (path_result, "_estimator_MDS.png", sep = ""),  1700, 1500)
	par( mar=c(8,8,8,8))
	plot (fit$points[,1], fit$points[,2], main="", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, pch = 20, col = color_ALSL)	
	#text (fit$points[,1], fit$points[,2], labels = rownames (data), col = as.integer(color_point),  cex = 2.5)
	abline(h=0,v=0)
	dev.off()

	svg (file =paste (path_result, "_estimator_MDS.svg", sep = ""), 25, 25)
	plot (fit$points[,1], fit$points[,2], main="", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, pch = 20, col = color_ALSL)	
	abline(h=0,v=0)
	dev.off()

	###############################
	# color by type druggabbility #
	###############################

	png (paste (path_result, "_drugg_MDS.png", sep = ""), 1700, 1500)
	par(mar=c(8,8,8,8))
	color_point[which (color_point == "#01B0F0")] = 2
	color_point[which (color_point == "#12A003")] = 1
	plot(fit$points[,1], fit$points[,2], main="", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, pch = 20, col = color_point)	
	abline(h=0,v=0)
	dev.off()

	svg (file =paste (path_result, "_drugg_MDS.svg", sep = ""), 25, 25)
	plot(fit$points[,1], fit$points[,2], main="", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, pch = 20, col = color_point)		
	abline(h=0,v=0)
	dev.off()	
	

}





################
#     MAIN     #
################

args <- commandArgs(TRUE)
path_dataset1 = args[1]
path_dataset2 = args[2]
path_correspondance = as.character(args[3])
path_result = args[4]

#print (path_dataset1)
#print (path_dataset2)

if (path_correspondance != "0"){
	matrix_correspondance = read.table (path_correspondance, header = FALSE, sep = "\t")
}else {
	matrix_correspondance = 0


}


data1_open = openData (path_dataset1, 0)
descriptor_data1 = data1_open[[2]]


data2_open = openData (path_dataset2, 0)
descriptor_data2 = data2_open[[2]]


descriptor_identic = c(intersect (descriptor_data1, descriptor_data2), "drugg")

data1 = subset(data1_open[[1]], select = descriptor_identic)
data2 =  subset(data2_open[[1]], select = descriptor_identic)


drug1 = data1$drugg + 1
drug2 = data2$drugg

data1 = data1[,-which(colnames (data1) == "drugg")]
data2 = data2[,-which(colnames (data2) == "drugg")]

#print (dim (data1))
#print (dim (data2))

#MDSCor (data1, data2, path_result, drug1, drug2, matrix_correspondance)
MDSDist (data1, data2, path_result, drug1, drug2, matrix_correspondance)




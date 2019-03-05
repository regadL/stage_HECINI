#!/usr/bin/env Rscript

# By BORREL Alexandre
# 04-2012

source ("tool.R")


MDS = function (path_file, col.desc){

	data = read.table (path_file, sep = "\t", header = TRUE)
	data = na.omit (data)
	data = data[,-which (colnames(data) == "drugg")]

	MC = cor(data)
	dist1 = 1-MC
	
	fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)

	png (paste (path_file, "_MDS.png", sep = ""), 2480, 1700)
	par( mar=c(12,12,12,12))
	fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
	c = as.vector (col.desc[,colnames (data)])
	plot (fit$points[,1], fit$points[,2], main="MDS Descriptor", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, type = "n")
	text (fit$points[,1], fit$points[,2], labels = colnames (data), col = as.integer(c),  cex = 2.5)
	legend("topleft", col=c(3,4,6), legend = c("desc RADI","desc perso","desc Fpocket"),pch=c(26,26,26),lty=c(1,1,1), cex = 1.5)  
	
	dev.off()
}



########
# MAIN #
########

args <- commandArgs(TRUE)
path_file = args[1]

d.col = colorACPByTypeOfDescriptors ()

MDS (path_file, d.col)


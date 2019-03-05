#!/usr/bin/env Rscript

# By BORREL Alexandre
# 04-2012

source ("tool.R")
require(graphics)



ACP = function (data, color_point, path_png, name, col.desc, type_plot){
	
	# normalization
	data = scale (data)

	# data correlation
	data.cor=cor(data)
	
	# eigen vector
	data.eigen=eigen(data.cor)

	lambda = data.eigen$values
	var_cap = cumsum(lambda)/sum(lambda)*100
	
	cp = data.eigen$vectors
	rownames (cp) = colnames (data)
	colnames (cp) = colnames (data)
	data_plot = as.matrix(data)%*%cp
	#print (data_plot)
	factor = factorACP (data_plot, cp) 
	if (type_plot == "ps"){
		postscript (file = paste(path_png, ".ps", sep = ""), width = 0, height = 0)
	}
	else{
		png (paste(path_png, ".png", sep = ""), 1000, 1000)
        par(mar=c(8,8,8,8))
	}
	par(mar=c(8,8,1,8))
	plot(data_plot[,1],data_plot[,2], col = color_point, pch=19, xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 1.8, cex.main = 2, cex = 2)
	abline(h=0,v=0)

 	#legend("topleft", col=c(2,1,3,4,6), legend = c("Drugg", "NoDrugg","desc RADI","desc perso","desc Fpocket"),pch=c(19,19,26,26,26),lty=c(0,0,1,1,1), cex = 1.5)  
	warnings ()	
	dev.off()

	if (type_plot == "ps"){
		postscript (file = paste(path_png, "_descriptor.ps", sep = ""), width = 0, height = 0)
	}else{
		png (paste(path_png, "_descriptor.png", sep = ""), 1000, 1000)
        par(mar=c(8,8,8,8))
	}

	par(mar=c(8,8,1,8))
	plot(data_plot[,1],data_plot[,2], col = color_point, pch=19, xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 1.8, cex.main = 2, type = "n") #, xlim = c(-10,10), ylim = c(-8,8))
	color_arrow = as.vector (col.desc[rownames(cp)])

	arrows(0,0,cp[,1]*factor,cp[,2]*factor, col= as.integer(color_arrow))
	text(cp[,1]*factor,cp[,2]*factor, rownames(cp),  col= as.integer(color_arrow), cex = 1.5)
	abline(h=0,v=0)
 	#legend("topleft", col=c(2,1,3,4,6), legend = c("Drugg", "NoDrugg","desc RADI","desc perso","desc Fpocket"),pch=c(19,19,26,26,26),lty=c(0,0,1,1,1), cex = 1.5)  
	warnings ()	
	dev.off()




	png (paste(path_png, "_PDB.png", sep = ""), 1000, 1000)
	par(mar=c(8,8,1,8))
	plot(data_plot[,1],data_plot[,2], col = color_point, main = "ACP dataset PDB", pch=20, xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
	text (data_plot[,1],data_plot[,2], col = color_point, label = rownames (data), cex = 1.6)
	dev.off()



}



################
#     MAIN     #
################

args <- commandArgs(TRUE)
path_dataset = args[1]
name = as.character (args[2])


d.col = colorACPByTypeOfDescriptors ()

data = openData (path_dataset, 0)[[1]]

if (is.integer0(which(colnames(data) == "drugg")) == FALSE){
	col_point = data$drugg + 1
	data = data[,-which(colnames (data) == "drugg")]
}else{
	col_point = rep(1,dim(data)[2])
}

d.col = CheckColorVector (colnames (data), d.col)

#ACP(data, col_point, paste (path_dataset, "_ACP", sep = ""), name, d.col, "ps")
ACP(data, col_point, paste (path_dataset, "_ACP", sep = ""), name, d.col, "png")

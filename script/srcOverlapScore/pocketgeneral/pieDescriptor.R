#!/usr/bin/env Rscript

# By BORREL Alexandre
# 10-2013


source ("tool.R")
library (MASS)




###########
## MAIN  ##
###########

args <- commandArgs(TRUE)
path_filin = args[1]


d.col = colorACPByTypeOfDescriptors()
d_g = read.table (path_filin, row.names = 1)

#print (dim(d_g))


if (dim(d_g)[2] != 1){

	png (paste(path_filin, dim(d_g)[2], "_pie.png", sep = ""),1000*dim(d_g)[2], 600)
	par (mfrow = c(1,dim(d_g)[2]))

	for (i in seq (dim(d_g)[2])){
		d = d_g[- which(d_g[,i] == 0),]
		d = orderByType (d.col, d)
		cpie = d.col[rownames(d)]
		pie(d[,i],labels = rownames (d), col = as.double(cpie)-1)


	}

	dev.off()

}else {


	d = d_g 
	d = orderByTypeOneCol (d.col, d)
	cpie = d.col[,names(d)]

	png (paste(path_filin, "_pie.png", sep = ""), 800, 500)
	pie(d,labels = names (d), col = as.double(cpie)-1)
	dev.off()



}




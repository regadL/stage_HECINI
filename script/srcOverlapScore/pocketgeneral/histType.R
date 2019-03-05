#!/usr/bin/env Rscript

# By BORREL Alexandre
# 05-2012


histNameProtein = function (data_plot, path_file){
	
	
	data_plot = data_plot[order(data_plot[,2]),]

	if (dim(data_plot)[1]> 200){
		print (as.integer((dim (data_plot)[1] / 200 )))
		png (paste (path_file, ".png", sep = ""), as.integer((dim (data_plot)[1] / 200 )* 3508), 2480)
	}
	else{
		png (paste (path_file, ".png", sep = ""), 3508, 2480)	
	}
	par( mar=c(75,10,5.5,1.5))
	color = data_plot[,3]
	for (i in seq (1, length (color))){
		if (color[i] == 2){
			color = color - 1
			break
		}
	}
	
	barplot(data_plot[,2], names.arg = data_plot[,4], las = 2, ylim = c(0,1), cex.names= 2.8,cex.axis= 3, col = color, cex.lab = 2)
	grid(NA, 10, col = 1, lwd = 2) 

	dev.off ()

}

barplot2Col = function (data, path_file_data){

	png (paste (path_file_data, ".png", sep = ""), 3000, 1500)
	par( mar=c(20,10,10,10))
	data=data[order(data[,2],decreasing = T),]

	data = data[1:200,]

	barplot (data[,2], names.arg = data[,1], las = 2)

	dev.off()
}



############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file = args[1]
data = read.table (path_file, sep = "\t", header = FALSE)

if (dim(data)[2] == 2){
	barplot2Col (data, path_file)
}else{
	histNameProtein (data, path_file)
}

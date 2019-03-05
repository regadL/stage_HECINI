#!/usr/bin/env Rscript

# By BORREL Alexandre
# 05-2012

source ("tool.R")


ACP = function (data1, data2, path_result, color_point, color_point2, matrix_correspondance){
	#print (data)	

	color_point2[which(color_point2 == 1)] = "#01B0F0"
	color_point2[which(color_point2 == 0)] = "#12A003"

	#Scale
	data1 = scale (data1)
	data2 = scale (data2)

	#Coord data 1
	data1.cor=cor(data1)
	data1.eigen=eigen(data1.cor)
	lambda1 = data1.eigen$values
	var_cap = cumsum(lambda1)/sum(lambda1)*100
	cp1 = data1.eigen$vectors
	rownames (cp1) = colnames (data1)
	colnames (cp1) = colnames (data1)
	data_plot1 = as.matrix(data1)%*%cp1

	# Coord data 2
	data2.cor=cor(data2)
	data2.eigen=eigen(data2.cor)
	

	lambda2 = data2.eigen$values
	cp2 = data2.eigen$vectors
	rownames (cp2) = colnames (data2)
	colnames (cp2) = colnames (data2)
	data_plot2 = as.matrix(data2)%*%cp2


	#png (paste (path_result, ".png", sep = ""), 2480, 3508)
	#par(mar=c(8,8,8,8))
	#par (mfrow = c(2,1))
	#plot(data_plot1[,1],data_plot1[,2], col = color_point, main = "ACP dataset", pch=20, xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
	#abline(h=0,v=0)
	#points (data_plot2[,1],data_plot2[,2], col =color_point2, pch=17, cex = 4)
 	#legend("topleft", col=c(2,1,"#F100F1","#12A003"), legend = c("Drug data1", "NoDrug data1","Drug data2", "NoDrug data2") ,pch=c(20,20,17,17),lty=c(0,0,0,0), cex = 3.5)  
	
	#plot(data_plot1[,1],data_plot1[,2], col = color_point, main = "ACP dataset", type = "n", xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
	#abline(h=0,v=0)
	#text(data_plot1[,1],data_plot1[,2], col = color_point, labels = rownames (data_plot1), cex = 2.5)
	#points (data_plot2[,1],data_plot2[,2],  type = "n")
	#text(data_plot2[,1],data_plot2[,2], col = color_point2, labels = rownames (data_plot2), cex = 2.5)
 	#legend("topleft", col=c(2,1,"#F100F1","#12A003"), legend = c("Drug data1", "NoDrug data1","Drug data2", "NoDrug data2") ,pch=c(20,20,20,20),lty=c(0,0,0,0), cex = 3.5) 

	#warnings ()	
	#dev.off()


	data = rbind (data1,data2)

	#print (dim (data1))
	#print (dim (data2))
	#print (dim (data))

	color_point = append (color_point, color_point2)
	#print (length (color_point))
	#shape = rep (20, length (color_point))
	#shape = append (shape, rep (17, length (color_point2)))	
	
	#data = scale (data)

	data.cor=cor(data)
	data.eigen=eigen(data.cor)
	lambda = data.eigen$values
	var_cap = lambda/sum(lambda)*100
	cp = data.eigen$vectors
	rownames (cp) = colnames (data)
	colnames (cp) = colnames (data)
	data_plot = as.matrix(data)%*%cp

	png (paste (path_result, "_PDB.png", sep = ""), 1700, 1500)
	factor = factorACP (data_plot, cp)
	#print (factor)
	col.desc = colorACPByTypeOfDescriptors ()
	#print (col.desc)
	#print (col.desc)	
	#print (rownames(cp))
	color_arrow = col.desc[rownames(cp)]
	#print (rapply(color_arrow, cbind))
	#print (color_arrow)
	par(mar=c(8,8,8,8))
	#plot(data_plot[,1],data_plot[,2], type = "n")
	plot(data_plot[,1],data_plot[,2], col = color_point, pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = "CP1", ylab = "CP2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
	text (data_plot[,1],data_plot[,2], col = color_point, label = rownames (data), cex = 1.6)	
	#points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)	
	abline(h=0,v=0)
	#arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = as.integer(color_arrow), lwd = 3 )
	#text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = as.integer(color_arrow), cex = 2.5)
 	#legend("topleft", col=c(2,1,"#01B0F0","#12A003"), legend = c("Drug data1", "NoDrug data1","Drug data2", "NoDrug data2") ,pch=c(20,20,20,20),lty=c(0,0,0,0), cex = 3.5)  
	warnings ()	
	dev.off()


	png (paste (path_result, "_color.png", sep = ""), 1700, 1500)
	factor = factorACP (data_plot, cp)
	color_ALSL = color_point

	color_ALSL[which(color_ALSL == 1)] = 3
	color_ALSL[which(color_ALSL == 2)] = 3
	color_ALSL[which(color_ALSL == "#01B0F0")] = 4
	color_ALSL[which(color_ALSL == "#12A003")] = 4	
	#print (color_ALSL)
	#print (factor)
	col.desc = colorACPByTypeOfDescriptors ()
	#print (col.desc)	
	#print (rownames(cp))
	color_arrow =col.desc[rownames(cp)]
	#print (color_arrow)
	par(mar=c(8,8,8,8))
	#plot(data_plot[,1],data_plot[,2], type = "n")
	plot(data_plot[,1],data_plot[,2], col = color_ALSL, pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = "CP1", ylab = "CP2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
	#text (data_plot[,1],data_plot[,2], col = color_point, label = rownames (data), cex = 1.6)	
	#points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)	
	abline(h=0,v=0)
	#arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = as.integer(color_arrow), lwd = 3 )
	#text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = as.integer(color_arrow), cex = 2.5)
 	#legend("topleft", col=c(2,1,"#01B0F0","#12A003"), legend = c("Drug data1", "NoDrug data1","Drug data2", "NoDrug data2") ,pch=c(20,20,20,20),lty=c(0,0,0,0), cex = 3.5)  
	warnings ()	
	dev.off()


	svg (paste (path_result, "_color.svg", sep = ""), 25, 25)
	factor = factorACP (data_plot, cp)
	color_ALSL = color_point

	color_ALSL[which(color_ALSL == 1)] = 3
	color_ALSL[which(color_ALSL == 2)] = 3
	color_ALSL[which(color_ALSL == "#01B0F0")] = 4
	color_ALSL[which(color_ALSL == "#12A003")] = 4	
	#print (color_ALSL)
	#print (factor)
	col.desc = colorACPByTypeOfDescriptors ()
	#print (col.desc)	
	#print (rownames(cp))
	color_arrow = col.desc[rownames(cp)]
	#print (color_arrow)
	par(mar=c(8,8,8,8))
	#plot(data_plot[,1],data_plot[,2], type = "n")
	plot(data_plot[,1],data_plot[,2], col = color_ALSL, pch=20, main = "", xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 6, cex.main = 4, cex.axis = 1.75, cex = 6)
	#text (data_plot[,1],data_plot[,2], col = color_point, label = rownames (data), cex = 1.6)	
	#points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)	
	abline(h=0,v=0)
	#arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = as.integer(color_arrow), lwd = 3 )
	#text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = as.integer(color_arrow), cex = 2.5)
 	#legend("topleft", col=c(2,1,"#01B0F0","#12A003"), legend = c("Drug data1", "NoDrug data1","Drug data2", "NoDrug data2") ,pch=c(20,20,20,20),lty=c(0,0,0,0), cex = 3.5)  
	warnings ()	
	dev.off()




	png (paste (path_result, "_descriptor.png", sep = ""), 1700, 1500)
	par(mar=c(8,8,8,8))
	plot(data_plot[,1],data_plot[,2], col = color_point, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = "CP1", ylab = "CP2", pch=20, cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type = "n")
	#points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)	
	abline(h=0,v=0)
	arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = as.character(color_arrow), lwd = 3 )
	text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = as.character(color_arrow), cex = 2.5)
	dev.off()


	svg (file = paste (path_result, "_descriptor.svg", sep = ""), 25, 25)
	par(mar=c(8,8,8,8))
	plot(data_plot[,1],data_plot[,2], col = color_point, main = "", xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 6, cex.main = 4, cex.axis = 1.75, cex = 6, type = "n")
	#points (data_plot[,1][length (color_point):dim(data_plot)[1]],data_plot[,2][length (color_point):dim(data_plot)[1]], pch=17, cex = 4, col = color_point2)	
	abline(h=0,v=0)
	arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = as.character(color_arrow), lwd = 3 )
	text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = as.character(color_arrow), cex = 3.5)
	dev.off()	


	png (paste (path_result, "_points.png", sep = ""), 1700, 1500)
	par(mar=c(8,8,8,8))
	color_point[which (color_point == "#01B0F0")] = 2
	color_point[which (color_point == "#12A003")] = 1
	plot(data_plot[,1],data_plot[,2], col = color_point, pch=c(rep(20,115), rep (8, 112)), main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = "CP1", ylab = "CP2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
	abline(h=0,v=0)
	dev.off()

	svg (file = paste (path_result, "_points.svg", sep = ""), 25, 25)
	par(mar=c(8,8,8,8))
	color_point[which (color_point == "#01B0F0")] = 2
	color_point[which (color_point == "#12A003")] = 1
	plot(data_plot[,1],data_plot[,2], col = color_point, pch=20, main = "", xlab = paste("CP1: ", signif (var_cap[1], 4), "%", sep = ""), ylab = paste("CP2: ", signif (var_cap[2], 4), "%", sep = ""), cex.lab = 6, cex.main = 4, cex.axis = 1.75, cex = 6)
	abline(h=0,v=0)
	dev.off()

	color_drug = color_point
	
	if (matrix_correspondance != "0"){	
		png (paste (path_result, "color_family.png", sep = ""), 1700, 1500)
		par(mar=c(8,8,8,8))

		matrix_correspondance = cbind(matrix_correspondance, rainbow(dim(matrix_correspondance)[1]))
		PDB = append (as.character(matrix_correspondance[,1]),as.character(matrix_correspondance[,2]))
		color_family = append (rainbow(dim(matrix_correspondance)[1]),rainbow(dim(matrix_correspondance)[1]))
		color_point = cbind (PDB, color_family)		
		rownames(color_point) = color_point[,1]
		
		color_point = color_point[rownames(data_plot),]
	
		plot(data_plot[,1],data_plot[,2], col = color_point[,2], pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = "CP1", ylab = "CP2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4)
		abline(h=0,v=0)
		dev.off()


		svg (file = paste (path_result, "_fig_article.svg", sep = ""),25,25)
		par(mar=c(8,8,8,8))
		plot(data_plot[,1],data_plot[,2], col = color_point[,2],pch=20, main = paste (var_cap[1],var_cap[2], sep = "_" ), xlab = "CP1", ylab = "CP2", cex.lab = 4, cex.main = 4, cex.axis = 1.75, cex = 4, type="n")
		#print (rownames(data_plot))	
		text(data_plot[,1],data_plot[,2], labels = c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","V", "VII", "X", "III", "I", "XIV", "XIII", "VI", "IX", "XI", "IV", "XII", "VIII", "II", "XV", "XVII", "XVI"), col = color_drug, cex=2.5)
		#arrows (0,0,cp[,1]*factor,cp[,2]*factor, col = as.character(color_arrow), lwd = 3 )
		#text (cp[,1]*factor, cp[,2]*factor,rownames(cp), col = as.character(color_arrow), cex = 2.5)
		abline(h=0,v=0)
		dev.off()
	}
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

#name_2 = rownames(data2)
#print (name_2) 
#del_element = match( name_2,rownames (data1))


#data1 = data1[-del_element,]

drug1 = data1$drugg + 1
drug2 = data2$drugg

data1 = data1[,-which(colnames (data1) == "drugg")]
data2 = data2[,-which(colnames (data2) == "drugg")]

#print (dim (data1))
#print (dim (data2))

ACP (data1, data2, path_result, drug1, drug2, matrix_correspondance)





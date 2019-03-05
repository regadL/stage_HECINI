#!/usr/bin/env Rscript


# By BORREL Alexandre
# 10-2013 (redit for save model)

source ("tool.R")
require(plotrix)
library (MASS)


modelLDA = function (desc){
	
	Y_out = changeList (desc$drugg)
	desc_lda = desc[,-which(colnames(desc) == "drugg")] # remove drugg column

	res.lda = lda(x=desc_lda,grouping=Y_out) # model
	
	return (res.lda)
	
}


barplotDescriptor = function(res.lda, name_file, data_train){

	aa = dim (data_train)
	#print (attributes(res.lda))
	#print ("ggggg")
	#print (res.lda$scaling)
	#print ("dddddd")
	# coeficent LDA
	coef = res.lda$scaling[,1]
				
	# Standardized
	coef = normalizationScalingLDA (coef, data_train)

	vcol = rep(2,length=length(coef))
	names(vcol) = names(coef)
	vcol[which(coef<0)]=4
	coef.sort = sort (abs(coef),decreasing = TRUE)
	print (coef.sort)
	png (paste (name_file, ".png", sep = ""), 2000, 2000)
	par( mar=c(45,10,5,5))
	barplot(coef.sort, names.arg = names(coef.sort), las=2, cex.names = 2.75, col = vcol[names(coef.sort)], main = "", cex.main = 3, cex.axis = 3, ylab = "Significativité des descripteurs (%)", cex.lab = 3)
	#legend("topright",legend=c("coef >0","coef <0"), col=c(2,4),lty=2, cex = 2.75)
	dev.off ()

}






############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file_train = args[1]


# open for LDA
data_open = openData (path_file_train, 0)
descriptor = data_open[[2]]
data_train = data_open[[1]]


model.lda = modelLDA (data_train)

barplotDescriptor (model.lda, paste(path_file_train, "_desc", sep = ""), data_train)
save(model.lda, file = paste(path_file_train, ".Rdata", sep = ""))

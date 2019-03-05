#!/usr/bin/env Rscript


# By BORREL Alexandre
# 06-2012

source ("tool.R")
library (rpart)
library (MASS)


cartTree = function (data_descriptor, path_tree){

	model_cart = rpart (drugg ~., data = data_descriptor, control=rpart.control(minsplit=3), method = "class")

	png (paste (path_tree, ".png", sep = ""), 2000, 1700)
	plot (model_cart, cex.main = 2)
	text (model_cart, use.n=TRUE, cex = 1.75)
	dev.off()
}



#################################
#           MAIN                #
#################################

args <- commandArgs(TRUE)
path_file_global = args[1]
path_file_tree = args[2]


res.open = openData (path_file_global, 0)
data_global = res.open[[1]]


cartTree (data_global, path_file_tree)

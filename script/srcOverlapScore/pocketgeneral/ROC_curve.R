#!/usr/bin/env Rscript


# By BORREL Alexandre
# 09-2012

source ("tool.R")
library (MASS)
source ("predictive_power.R")

ldaROC = function (data_in, model_LDA){
	
	Y = as.factor(data_in$drugg)
	list_real = rep("d",length(Y))
	list_real[which(Y==0)]="nd"

	data_in = data_in[,-which(colnames(data_in) == "drugg")] # remove drugg column
	


	v_predict_proba = predict (model_LDA, data_in)$posterior # prediction


	vect_se = NULL
	vect_sp = NULL

	for (prob in seq(0.1,1,0.05)){
		vect_predict = generateVect (v_predict_proba, prob)
		rate_temp = calculTaux2 (vect_predict, list_real)

		vect_se = c(vect_se,sensibility(rate_temp[1], rate_temp[4]))
		vect_sp = c(vect_sp,1-specificity(rate_temp[2], rate_temp[3]))
	}
	return (list (vect_se, vect_sp))
}




generateLDA = function (data_train){
	
	Y = as.factor(data_train$drugg)
	Y2 = rep("d",length(Y))
	Y2[which(Y==0)]="nd"
	
	data_train = data_train[,-which(colnames(data_train) == "drugg")] # remove drugg column

	res.lda = lda(x=data_train,grouping=Y2) # LDA model
	return (res.lda)
}






############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file_train = args[1]
path_file_test = args[2]
path_png = args[3]

data_open = openData (path_file_train, 0)
desc = data_open[[1]]
descriptor = data_open[[2]]



desc_test = openData (path_file_test, 0)[[1]]

# model generation
model_LDA = generateLDA(desc)


vec_curve_train = ldaROC (desc, model_LDA)
vec_curve_test = ldaROC (desc_test, model_LDA)

png (path_png)

plot (vec_curve_train[[2]], vec_curve_train[[1]], pch = 19, col = "red", xlab = "1-Sp", ylab = "Se", xlim = c(0,1), ylim = c(0,1))
points(vec_curve_test[[2]], vec_curve_test[[1]], pch = 4, col = "blue")

graphics.off()








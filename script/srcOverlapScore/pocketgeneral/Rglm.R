#!/usr/bin/env Rscript


# By BORREL Alexandre
# 10-2013

library (MASS)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")



AIC = function (data_train, path_file_global){

	resAIC=NULL
	respval=NULL

	for (i in 1 : (dim(data_train)[2]-1)) {
		modele = glm(drugg~data_train[,i],data=data_train,family=binomial)
		resAIC = c(resAIC,modele$aic)   ###stocke la valeur de l'aic du modele
		respval = c(respval ,summary(modele)$coefficients[2,4])  ###stocke la pvalue du coefficient de la variable
	}
	names(resAIC) = colnames(data_train)[-dim(data_train)[2]]
	names(respval) = colnames(data_train)[-dim(data_train)[2]]


	resAIC.sort = sort(resAIC, decreasing = T)
	respval.sort = respval[names(resAIC.sort)]

	png (paste (path_file_global, "_residuel.png", sep = ""), 1600, 2000)
	par (mar = c(40,5,5,5))
	barplot(resAIC.sort, names.arg = names(resAIC.sort), las=2, cex.names = 2.75, main = "res AIC", cex.main = 3, cex.axis = 3, ylab = "Descriptors residuel", cex.lab = 3)
	dev.off()

	png (paste (path_file_global, "_coef.png", sep = ""), 1600, 2000)
	par (mar = c(40,5,5,5))
	barplot(respval.sort, names.arg = names(respval.sort), las=2, cex.names = 2.75, main = "Coef AIC", cex.main = 3, cex.axis = 3, ylab = "Descriptors coefficient", cex.lab = 3)
	dev.off ()

}



GLMLOO = function (data){

	# class
	
	# sample	
	i_sample = sample (dim(data)[1])
	data = data[i_sample,]

	print ("--------------------------------")
	print ("--------------LOO---------------")

		
	y_real = data[,which(colnames(data)=="drugg")]
	
	v_predict = NULL
	for (i_line in seq(dim(data)[1])){
		data_train = data[-i_line,]
		data_test = data[i_line,]
		model_glm = glm (drugg ~., data = data_train)
		y_predicted = predict (model_glm, newdata = data_test[,-which(colnames (data_test)=="drugg")], type="response")
		v_predict = append (v_predict,y_predicted)
	}
	# chane proba 0 vs 1
	v_predict[which(v_predict >= 0.5)] = 1
	v_predict[which(v_predict < 0.5)] = 0

	out = qualityPredict (list(changeList(v_predict)), list(changeList(y_real)))

}


GLMTrainTest = function (data_train, data_test){

	# sample
	i_sample = sample (dim(data_train)[1])
	data_train = data_train[i_sample,]

	i_sample = sample (dim(data_test)[1])
	data_test = data_test[i_sample,]	

	# real
	y_real_test = data_test[,which(colnames(data_test)=="drugg")]
	y_real_train = data_train[,which(colnames(data_train)=="drugg")]
	data_test = data_test[,-which(colnames(data_test)=="drugg")]

	fit = glm (drugg ~., data = data_train)
	
	# trainning dataset
	y_predicted_train = predict(fit, newdata = data_train, type="response")
	y_predicted_train[which(y_predicted_train >= 0.5)] = 1
	y_predicted_train[which(y_predicted_train < 0.5)] = 0

	print ("--------------------------------")
	print ("------------Data train----------")
	out = qualityPredict (list(changeList(y_predicted_train)), list(changeList(y_real_train)))

	# test dataset
	y_predicted_test = predict(fit, newdata = data_test, type="response")
	y_predicted_test[which(y_predicted_test >= 0.5)] = 1
	y_predicted_test[which(y_predicted_test < 0.5)] = 0

	print ("--------------------------------")
	print ("------------Data test----------")
	out = qualityPredict (list(changeList(y_predicted_test)), list(changeList(y_real_test)))	

}




#################################
#           MAIN                #
#################################

args <- commandArgs(TRUE)
path_file_global = args[1]
path_file_train = args[2]
path_file_test = args[3]
elimcor = args[4]


res.open = openData (path_file_global, elimcor)
data_global = res.open[[1]]
descriptor = res.open[[2]]
descriptor = append (descriptor, "drugg")

data_train = openData (path_file_train, elimcor)[[1]][,descriptor]

print ("Elimcor")
print (elimcor)

# AIC glm


# LOO for glm
GLMLOO (data_global)


if (path_file_train != 0){
	data_train = openData (path_file_train, elimcor)[[1]][,descriptor]
	data_test = openData (path_file_test, elimcor)[[1]][,descriptor]
	# AIC glm
	AIC (data_train, path_file_global)
	# model on train / test
	GLMTrainTest (data_train, data_test)
}


#!/usr/bin/env Rscript

# By BORREL Alexandre
# 05-2012


library (kernlab)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")
source ("learningMachine.R")



resultLOO = function (data, parameter){
	
	print ("----Leave one out----")
	y_real = data[,which(colnames(data)=="drugg")]
	y_predict = SVMLOO (data)
	
	out = qualityPredict (list(y_predict), list(changeList(y_real)))
}



SVMTEST = function (data_train,data_test, parameter){


	i_sample = sample (dim(data_train)[1])
	data_train = data_train[i_sample,]

	i_sample = sample (dim(data_test)[1])
	data_test = data_test[i_sample,]	

	y_real_test = data_test[,which(colnames(data_test)=="drugg")]
	y_real_train = data_train[,which(colnames(data_train)=="drugg")]
	data_test = data_test[,-which(colnames(data_test)=="drugg")]

	model_svm = ksvm(drugg~., data_train, kernel="rbfdot", kpar=list(sigma=parameter[1]),C=parameter[2], type = "C-svc")
	
	# trainning dataset
	y_predicted_train = predict(model_svm, data_train)
	print ("--------------------------------")
	print ("------------Data train----------")
	out = qualityPredict (list(changeList(y_predicted_train)), list(changeList(y_real_train)))

	# test dataset
	y_predicted_test = predict(model_svm, data_test)
	print ("--------------------------------")
	print ("------------Data test----------")
	out = qualityPredict (list(changeList(y_predicted_test)), list(changeList(y_real_test)))	
}


##############
##   MAIN   ##
##############

args <- commandArgs(TRUE)
path_file_train = args[1]
path_file_test = as.character(args[2])

open_data_train = openData (path_file_train, 0)
data_train = open_data_train[[1]]
data_train = data_train [sample (dim(data_train)[1]),]

print (colnames (data_train))

if (path_file_test != 0){
	open_data_test = openData (path_file_test, 0)
	data_test = open_data_test[[1]]
	data_test = data_test [sample (dim (data_test)[1]),]
	data_global = rbind (data_train, data_test)
}else {
	data_global = data_train
}

resultLOO (data_global)

if (path_file_test != "0"){
	parameter = optimization (data_train)
	print ("-> PARAMETERS TRAIN<-")
	print (paste ("C ->", parameter[1], " ", "G ->", parameter[2], sep = ""))
	SVMTEST (data_train, data_test, parameter)
}

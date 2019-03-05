#!/usr/bin/env Rscript


# By BORREL Alexandre
# 10-2013



library (randomForest)
library (MASS)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")





parametersGrid = function (lntree, lmtry, data_train, path_global){
	i_sample = sample (dim(data_train)[1])
	data_train = data_train[i_sample,]
	y_real = data_train[,which(colnames(data_train)=="drugg")]	

	i = 0
	grid = data.frame ()
	for (ntree in lntree){
		i = i + 1
		j = 0
		for (mtry in lmtry){
			j = j + 1
			fit = randomForest ( as.factor (drugg)~., data = data_train, mtry=mtry, ntree = ntree , type = "class")
			l_predict = predict (fit, dataNew = data_train[,-which(colnames (data_train)=="drugg")], type = "class")
			
			# control NA in prediction
			if (is.integer0(which(is.na(as.character(l_predict)))) == FALSE){
				l_predict[which(is.na(l_predict))] = 1
			}
			
			# R conversion 
			rate = calculTaux2  (changeList(as.double(l_predict)-1),changeList(y_real))
			mcc = MCC (rate[1], rate[2], rate[3], rate[4])
			grid[i,j] = mcc
		}
	}
	colnames (grid) = lmtry
	rownames (grid) = lntree

	write.table (grid, paste(path_global, ".grid", sep = ""))

	return (list(rownames (grid)[which(grid==max(grid), arr.ind=T)[1]],colnames (grid)[which(grid==max(grid), arr.ind=T)[2]] ))
}


RFLOO = function (data, parameters){

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
		model_RF = randomForest (as.factor (drugg) ~., data = data_train, mtry = as.integer(parameters[[2]]), ntree = as.integer(parameters[[1]]), method = "class")
		y_predicted = predict (model_RF, newdata = data_test[,-which(colnames (data_test)=="drugg")], type="class")
		v_predict = append (v_predict,y_predicted)
	}
	# -1 predict because predict in 1vs2
	out = qualityPredict (list(changeList(v_predict-1)), list(changeList(y_real)))
}


RFTrainTest = function (data_train, data_test, parameters, path_file_global){

	# sample
	i_sample = sample (dim(data_train)[1])
	data_train = data_train[i_sample,]

	i_sample = sample (dim(data_test)[1])
	data_test = data_test[i_sample,]	

	# real
	y_real_test = data_test[,which(colnames(data_test)=="drugg")]
	y_real_train = data_train[,which(colnames(data_train)=="drugg")]
	data_test = data_test[,-which(colnames(data_test)=="drugg")]

	fit = randomForest (as.factor(drugg)~., data = data_train, mtry = as.integer(parameters[[2]]), ntree = as.integer(parameters[[1]]), type = "class", importance = TRUE)
	png(paste(path_file_global, ".png", sep = ""))
	plot (fit)
	dev.off()
	
	# trainning dataset
	y_predicted_train = predict(fit, newdata = data_train, type="class")
	print ("--------------------------------")
	print ("------------Data train----------")
	out = qualityPredict (list(changeList(y_predicted_train)), list(changeList(y_real_train)))

	# test dataset
	y_predicted_test = predict(fit, newdata = data_test, type="class")
	print ("--------------------------------")
	print ("------------Data test----------")
	out = qualityPredict (list(changeList(y_predicted_test)), list(changeList(y_real_test)))	

}

varSelect = function (data_train, parameter, path_file_global){

	MatVI = matrix(NA, ncol = 100, nrow = dim(data_train)[2])
	rownames(MatVI) = colnames(data_train)

	for (i in 1:100){
		Imp = randomForest(as.factor(drugg)~., data=data_train, ntree = as.integer(parameter[[1]]), mtry = as.integer(parameter[[2]]))$importance
		MatVI[rownames(Imp), i] = Imp[,1]	
	}

	###importance means
	VIMean = apply(MatVI, 1, mean)
	names(VIMean) = rownames(MatVI)
	VISD = apply(MatVI, 1, sd)
	names(VISD) = rownames(MatVI)

	VIMeanSort = sort(VIMean, decreasing =T)

	png(paste(path_file_global,"_importance_variableRF.png"), 800, 800)
	par( mar=c(10,4,4,4))
	plot(VIMeanSort, xaxt ="n", xlab="", pch = 19, ylab="Importance moyenne")
	axis(1,1:length(VIMeanSort), labels = names(VIMeanSort), las = 2, cex.axis = 0.6, cex = 2.75)
	dev.off()


	###select 6 variables
	return (names(VIMeanSort)[1:6])
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
print (colnames (data_global))

# grid parameters based on MCC
parameter = parametersGrid (c(10,50,100,200,500, 1000), c(1,2,3,4,5,10,15,20, 25, 30), data_train, path_file_global)


# if you want select variable -> first 6
desc_select = varSelect (data_train, parameter, path_file_global)
desc_select = append (desc_select, "drugg")

desc_select = colnames (data_train)

# parameter again with selected descriptors
#parameter = parametersGrid (c(10,50,100,200,500, 1000), c(1,2,3,4,5,10,15,20, 25, 30), data_train[,desc_select], path_file_global)
print (parameter)

RFLOO (data_global[,desc_select],parameter)


if (path_file_train != 0){
	data_train = openData (path_file_train, elimcor)[[1]][,descriptor]
	data_test = openData (path_file_test, elimcor)[[1]][,descriptor]
	RFTrainTest (data_train[,desc_select], data_test[,desc_select], parameter, path_file_global)
}


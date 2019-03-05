#!/usr/bin/env Rscript


# By BORREL Alexandre
# 11-2013

library (MASS)
library (plotrix)
library (e1071)



source ("predictive_power.R")
source ("tool.R")


generateModel = function(data_train, data_test, data_global, nb_cross){

	# real LOO
	i_sample_loo = sample (dim(data_global)[1])
	data_global = data_global[i_sample_loo,]
	class = as.factor(data_global$drugg)
	class_real_global = rep("d",length(class))
	class_real_global[which(class==0)]="nd"

	# real train	
	i_sample_train = sample (dim(data_train)[1])
	data_train = data_train[i_sample_train,]
	class = as.factor(data_train$drugg)
	class_real_train = rep("d",length(class))
	class_real_train[which(class==0)]="nd"

	# real test
	i_sample_test = sample (dim(data_test)[1])
	data_test = data_test[i_sample_test,]
	class = as.factor(data_test$drugg)
	class_real_test = rep("d",length(class))
	class_real_test[which(class==0)]="nd"

	
	print (data_train)
	
	# creat model ad optimization
	time_begin = proc.time ()
	#print( tune.svm (as.factor(drugg)~., data = data_train, gamma = 2^(-1:1), cost = 2^( 2:4), tunecontrol = tune.control (sampling = "fix"), cross = nb_cross ))
	model = tune.svm (drugg~., data = data_train, gamma = 2^(-1:1), cost = 2^( 2:4), tunecontrol = tune.control (sampling = "fix"), cross = nb_cross)$best.model
	time_optimization = proc.time () - time_begin 
	
	print ("##############")					
	print (time_optimization[[3]])
	print ("##############")


	# apply model
	# train
	v_predicted_train = predict(model, data_train)
	print (v_predicted_train)
	# test
	v_predicted_test=predict(model, data_test)
				
		
	# rate train
	rate_predict_train = calculTaux2 (v_predicted_train, class_real_train)	
	acc_train = accuracy(rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])
	se_train = sensibility(rate_predict_train[1], rate_predict_train[4])
	sp_train = sensibility(rate_predict_train[2], rate_predict_train[3])
					

	# rate loo

	v_predict = NULL
	for (i_line in seq(dim(data_global)[1])){
		data_train = data_global[-i_line,]
		data_test = data_global[i_line,]
		
		# creat model and 
		model_LOO =  tune.svm (as.factor(drugg)~., data = data_train, gamma = 2^(-1:1), cost = 2^( 2:4), tunecontrol = tune.control (sampling = "fix"), cross = nb_cross )$best.model
		y_predicted = predict (model_LOO, newdata = data_test, type="class")
		
		v_predict = append (v_predict,y_predicted)
		print (v_predict)
	}

	
	# rate loo
	rate_loo = calculTaux2 (v_predict, class_real_global)
	acc_loo = accuracy(rate_loo[1], rate_loo[2], rate_loo[3], rate_loo[4])
	se_loo = sensibility(rate_loo[1], rate_loo[4])
	sp_loo = sensibility(rate_loo[2], rate_loo[3])

	# rate test
	rate_predict_test = calculTaux2 (v_predicted_test, class_real_test)	
	acc_test = accuracy(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])	
	se_test = sensibility(rate_predict_test[1], rate_predict_test[4])
	sp_test = sensibility(rate_predict_test[2], rate_predict_test[3])

	# retrieve condition
	#if (acc_loo > 0.80 && acc_train > 0.80 && acc_test > 0.80){
	print ("descriptor")
	print (colnames (data_global))
	print ("Acc_loo --- Acc_train --- Acc_test")
	print (paste (acc_loo, acc_train, acc_test, sep = "---"))
	print ("Se_loo --- Sp_loo --- Se_train --- Sp_train --- Se_test --- Sp_test")
	print (paste (se_loo, sp_loo, se_train, sp_train, se_test, sp_test, sep = "---"))
	print ("**********************************************************************")
}





###########
## MAIN  ##
###########

args <- commandArgs(TRUE)
path_file_train = args[1]
path_file_test = args[2]
nb_cross_validation = as.integer(args[3])

##############
# Open Files #
##############

data_train_open = openData (path_file_train, 0)
data_train = data_train_open[[1]]
descriptor = data_train_open[[2]]
data_test = openData (path_file_test,0)[[1]]
data_global = rbind (data_train, data_test)


print ("For control")
print(dim (data_global))


########
# RUN  #
########




generateModel (data_train, data_test, data_global, nb_cross_validation)






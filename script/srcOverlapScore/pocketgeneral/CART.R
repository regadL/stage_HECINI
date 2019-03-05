#!/usr/bin/env Rscript


# By BORREL Alexandre
# 04-2012

library (rpart)
library (MASS)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")



cartTreeSelect = function (data_train, data_test, data_global, path_file){

	i_sample = sample (dim(data_train)[1])
	data_train = data_train[i_sample,]

	i_sample = sample (dim(data_test)[1])
	data_test = data_test[i_sample,]	

	y_real_test = data_test[,which(colnames(data_test)=="drugg")]
	y_real_train = data_train[,which(colnames(data_train)=="drugg")]
	data_test = data_test[,-which(colnames(data_test)=="drugg")]

	model_cart = rpart (drugg ~., data = data_global, method = "class")
	model_cart_control = rpart (drugg ~., data = data_global, method = "class", control=rpart.control(minsplit=3))
	
	# plot tree -> tree without control	
	png (paste (path_file, "_without_contol_CART.png", sep = ""), 2000, 1700)
	plot (model_cart_control, cex.main = 2)
	text (model_cart_control, use.n=TRUE, cex = 1.75)
	dev.off()

	png (paste (path_file, "_complexity.png", sep = ""))
	plotcp(model_cart_control)
	dev.off()


	# plot tree -> tree with control
	png (paste (path_file, "_with_contol_CART.png", sep = ""), 2000, 1700)
	plot (model_cart, cex.main = 2)
	text (model_cart, use.n=TRUE, cex = 1.75)
	dev.off()

}


cart = function (data_train, data_test, path_file){

	i_sample = sample (dim(data_train)[1])
	data_train = data_train[i_sample,]

	i_sample = sample (dim(data_test)[1])
	data_test = data_test[i_sample,]	

	y_real_test = data_test[,which(colnames(data_test)=="drugg")]
	y_real_train = data_train[,which(colnames(data_train)=="drugg")]
	data_test = data_test[,-which(colnames(data_test)=="drugg")]

	#model_cart = rpart (drugg ~., data = data_train, method = "class")
	model_cart = rpart (drugg ~., data = data_train, method = "class")
	
	# trainning dataset
	y_predicted_train = predict(model_cart, newdata = data_train, type="class")
	print ("--------------------------------")
	print ("------------Data train----------")
	out = qualityPredict (list(changeList(y_predicted_train)), list(changeList(y_real_train)))
	png (paste (path_file, "_train_CART.png", sep = ""), 2000, 1700)
	plot (model_cart, cex.main = 2)
	text (model_cart, use.n=TRUE, cex = 1.75)
	dev.off()

	# test dataset
	y_predicted_test = predict(model_cart, newdata = data_test, type="class")
	print ("--------------------------------")
	print ("------------Data test----------")
	out = qualityPredict (list(changeList(y_predicted_test)), list(changeList(y_real_test)))	

}


cartLOO = function (data, path_file){

	y_real = changeList (data$drugg)
	res.cart = rpart (drugg ~., data = data_global, method = "class") # for plot

	# plot all pocket
	png (paste (path_file, "_all_CART.png", sep = ""), 2000, 1700)
	plot (res.cart, cex.main = 2)
	text (res.cart, use.n=TRUE, cex = 1.75)
	dev.off()

	print ("--------------------------------")
	print ("--------------LOO---------------")

	i_sample = sample (dim(data)[1])
	data = data[i_sample,]	
	
	data_predic = data[,- which(colnames(data)=="drugg")]
	y_real = data[,which(colnames(data)=="drugg")]
	
	v_predict = NULL
	for (i_line in seq(dim(data)[1])){
		data_train = data[-i_line,]
		data_test = data_predic[i_line,]
		model_cart = rpart (drugg ~., data = data_train, control=rpart.control(minsplit=3), method = "class")
		y_predicted = predict (model_cart, newdata = data_test, type="class")
		v_predict = append (v_predict,y_predicted)
	}
	out = qualityPredict (list(changeList(v_predict)), list(changeList(y_real)))
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

print ("Elimcor")
print (elimcor)
print (dim(data_global[2]))
print (colnames (data_global))


cartLOO (data_global,path_file_global)


if (path_file_train != 0){
	data_train = openData (path_file_train, elimcor)[[1]][,descriptor]
	data_test = openData (path_file_test, elimcor)[[1]][,descriptor]

	cartTreeSelect (data_train, data_test, data_global, path_file_global)
	cart (data_train, data_test, path_file_global)
}

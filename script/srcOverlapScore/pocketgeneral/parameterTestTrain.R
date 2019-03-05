#!/usr/bin/env Rscript


# By BORREL Alexandre
# 04-2012

library (MASS)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")



modelLDA = function (data_train, data_test){


	
	# Y train
	Y = as.factor(data_train$drugg)
	Y2 = rep("d",length(Y))
	Y2[which(Y==0)]="nd"
	
	# Y test
	Y_test = as.factor(data_train$drugg)
	Y2_test = rep("d",length(Y_test))
	Y2_test[which(Y_test==0)]="nd"
	

	data_train = data_train[,-which(colnames(data_train) == "drugg")] # remove drugg column
	res.lda = lda(x=data_train,grouping=Y2) # LDA model
	v_predict_train = predict (res.lda, data_train) # prediction

	rate_predict_train = calculTaux2 (v_predict_train$class, Y2)
	acc_train = accuracy(rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])

	v_predict_test = predict (res.lda, data_test) 
	rate_predict_test = calculTaux2 (v_predict_test$class, Y2_test)
	acc_test = accuracy(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])

	return (c(acc_train, acc_test))
}




############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file_train = args[1]
path_file_test = args[2]
begin_elimcor = as.double(args[3])
step_elimcor = as.double(args[4])

print (path_file_train)
print (path_file_test)


acc_train = NULL
list_elimcor = NULL
acc_test = NULL
for (elimcor in seq (begin_elimcor, 1, step_elimcor)){

	list_elimcor = append (elimcor, list_elimcor)

	data_train_open = openData (path_file_train, elimcor)
	data_train = data_train_open[[1]]
	descriptor = data_train_open[[2]]
	data_test = openData (path_file_test,0)[[1]]
	data_test = data_test [,descriptor]

	acc_list = modelLDA (data_train, data_test)
	acc_train = append (acc_train, acc_list[1])
	acc_test = append (acc_test, acc_list[2])
}

png (paste(path_file_train, "_elimcor.png"), 2480, 1700)
plot (list_elimcor, acc_train, type = "l", col = "red")
par (new = TRUE)
plot (list_elimcor, acc_test, type = "l", col = "blue")
legend ("topleft", legend = c("Accuracy train","Accuracy test"),col = c("red", "blue" ), cex = 2.75, bty="y", lty=1)
dev.off()





#!/usr/bin/env Rscript

# By BORREL Alexandre
# 07-2012

source ("predictive_power.R")
library (MASS)
library (kernlab)
source ("tool.R")

################
# LDA learning #
################

# LDA -> give accuracy
LDAACC = function (data_train, data_test, descriptor_class){
	

	real_train = data_train[,which(colnames(data_train) == descriptor_class)]
	real_class_train = rep("d",length(real_train))
	real_class_train[which(real_train==0)]="nd"


	real_test = data_test[,which(colnames(data_test) == descriptor_class)]
	real_class_test = rep("d",length(real_test))
	real_class_test[which(real_test==0)]="nd"	
	

	data_train = data_train[,-which(colnames(data_train) == descriptor_class)] # remove drugg column
	# a verif
	data_test = data_test[,-which(colnames(data_test) == descriptor_class)] # remove drugg column
	
	# LDA run 
	res.lda = lda(x=data_train,grouping=real_class_train) # LDA model
	v_predict_test = predict (res.lda, data_test) # prediction on test
	v_predict_train = predict (res.lda, data_train)

	rate_train = calculTaux2 (v_predict_train$class, real_class_train)
	rate_test = calculTaux2 (v_predict_test$class, real_class_test)
	
	acc_train = accuracy(rate_train[1], rate_train[2], rate_train[3], rate_train[4])
	acc_test = accuracy(rate_test[1], rate_test[2], rate_test[3], rate_test[4])
	return (c(acc_train, acc_test))
}


LDALOO = function (data_train, descriptor_class){

	real_train = data_train[,which(colnames(data_train) == descriptor_class)]
	real_class_train = rep("d",length(real_train))
	real_class_train[which(real_train==0)]="nd"

	data_train = data_train[,-which(colnames(data_train) == descriptor_class)] # remove drugg column
	res.lda = lda(x=data_train,grouping=real_class_train, CV = TRUE) # leave-one-out

	rate_train = calculTaux2 (res.lda$class, real_class_train)
	acc_LOO = accuracy(rate_train[1], rate_train[2], rate_train[3], rate_train[4])

	return (acc_LOO)
}


ACPPredict = function (data, bad_predict, drug_vector, path_png, name){
	
	# ACP all points	
	color_point = rep (4, dim (data)[1])
	no_drug = which (drug_vector == "nd")
	no_drug_bad_predict = intersect (as.vector(no_drug), as.vector(bad_predict))
	color_point [bad_predict] = 2
	color_point[no_drug_bad_predict] = "#CC06D3"

	data.scale = scale (data)

	data.cor=cor(data.scale)
	data.eigen=eigen(data.cor)
	
	lambda = data.eigen$values
	var_cap = cumsum(lambda)/sum(lambda)*100
	
	cp = data.eigen$vectors
	rownames (cp) = colnames (data.scale)
	colnames (cp) = colnames (data.scale)
	data_plot = as.matrix(data.scale)%*%cp
	#print (data_plot)
	factor = factorACP (data_plot, cp)
	d.col = colorACPByTypeOfDescriptors ()
	c = as.vector (d.col[,rownames(cp)])

	png (paste(path_png, "_ACPGlobal.png", sep = ""), 2480, 3508)
	par(mfrow = c(2,1))	
	par(mar=c(8,8,8,8))

	# With druggable and non-druggable
	plot(data_plot[,1],data_plot[,2], type = "n", main = name, xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75)
	text (data_plot[,1],data_plot[,2], col = color_point, labels = drug_vector, cex = 2.2)	

	arrows(0, 0, cp[,1]*factor, cp[,2]*factor, col = as.integer(c))
	text(cp[,1]*factor, cp[,2]*factor, colnames(data.scale), col = as.integer(c), cex = 2.5)
	
	abline(h=0,v=0)
 	legend("topleft", col=c(4,2,"#CC06D3",3,4,6), legend = c("Good Predict", "Bad Druggable", "Bad Non-Druggable", "desc RADII3", "desc perso", "desc Fpocket"),pch=c(20,20,20,26,26,26),lty=c(0,0,0,1,1,1), cex = 3.5)
	
	# With complexes names
	plot(data_plot[,1],data_plot[,2], type = "n", main = name, xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75)
	text (data_plot[,1],data_plot[,2], col = color_point, labels = rownames (data.scale), cex = 2.2)	

	arrows(0, 0, cp[,1]*factor, cp[,2]*factor, col = as.integer(c))
	text(cp[,1]*factor, cp[,2]*factor, colnames(data.scale), col = as.integer(c), cex = 2.5)
	
	abline(h=0,v=0)
 	#legend("topleft", col=c(4,2,"#CC06D3",3,4,6), legend = c("Good Predict", "Bad Druggable", "Bad Non-Druggable", "desc RADII3", "desc perso", "desc Fpocket"),pch=c(20,20,20,26,26,26),lty=c(0,0,0,1,1,1), cex = 3.5)
	dev.off()

	#-------------------#
	# Only bad predicted#
	#-------------------#

	png (paste(path_png, "_ACPBad.png", sep = ""), 2480, 3508)
	par(mfrow = c(2,1))	
	par(mar=c(8,8,8,8))

	# ACP bad prediction
	data_bp = data[bad_predict,]
	data_bp = scale (data_bp)
	data_bp.cor=cor(data_bp)

	# case with cor nul (missing value or value = 0)	
	
	result <- try(data_bp.eigen<-eigen(data_bp.cor))
	if (result[1] == "Error in eigen(data_bp.cor) : infinite or missing values in 'x'\n"){
		return	()
	}

	data_bp.eigen=eigen(data_bp.cor)

	lambda = data_bp.eigen$values
	var_cap = cumsum(lambda)/sum(lambda)*100
	
	cp = data_bp.eigen$vectors
	rownames (cp) = colnames (data_bp)
	colnames (cp) = colnames (data_bp)
	data_plot = as.matrix(data_bp)%*%cp
	factor = factorACP (data_plot, cp)
	c = as.vector (d.col[,rownames(cp)])

	# With druggable and non-druggable
	plot(data_plot[,1],data_plot[,2], type = "n", main = paste (name, " bad prediction", sep = ""), xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75)
	text (data_plot[,1],data_plot[,2], col = 2, labels = drug_vector[bad_predict], cex = 2.2)	

	arrows(0, 0, cp[,1]*factor, cp[,2]*factor, col = as.integer(c))
	text(cp[,1]*factor, cp[,2]*factor, colnames(data_bp), col = as.integer(c), cex = 2.5)
	
	abline(h=0,v=0)
 	legend("topleft", col=c(2,3,4,6), legend = c("Bad Predict","RADI3","Personal","Fpocket"),pch=c(20,26,26,26),lty=c(0,1,1,1), cex = 3.5)

	# With complexes names
	plot(data_plot[,1],data_plot[,2], type = "n", main = paste (name, " bad prediction", sep = ""), xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75)
	text (data_plot[,1],data_plot[,2], col = 2, labels = rownames (data_bp), cex = 2.2)	

	arrows(0, 0, cp[,1]*factor, cp[,2]*factor, col = as.integer(c))
	text(cp[,1]*factor, cp[,2]*factor, colnames(data_bp), col = as.integer(c), cex = 2.5)
	
	abline(h=0,v=0)
 	legend("topleft", col=c(2,3,4,6), legend = c("Bad Predict","RADI3","Personal","Fpocket"),pch=c(20,26,26,26),lty=c(0,1,1,1), cex = 3.5)
	dev.off()
}


###################
#       SVM       #
###################


optimization = function (data_global){
	
	#i_sample = sample (dim(data_global)[1])
	#data_global = data_global[i_sample,]

	val_sigma = c(10^-4, 10^-3, 10^-2, 10^-1, 1, 10^1, 10^2, 10^3, 10^4)
	val_C = c(10^-4, 10^-3, 10^-2, 10^-1, 1, 10^1, 10^2, 10^3, 10^4)
	
	list_parameter = NULL
	for (sigma in 1:9){
		for (C in 1:9){
			y_real = data_global[,which(colnames(data_global)=="drugg")]
			data_predic = data_global[,- which(colnames(data_global)=="drugg")]
			
			model_svm = ksvm(drugg~., data_global, kernel="rbfdot", kpar=list(sigma=val_sigma[sigma]),C=val_C[C], type = "C-svc")
			y_predicted=predict(model_svm, data_predic)
			
			rate_svm = calculTaux2 (y_predicted, y_real)
			
			acc_temp = accuracy(rate_svm[1], rate_svm[2], rate_svm[3], rate_svm[4])
			list_parameter = rbind (list_parameter, c(val_sigma[sigma], val_C[C],acc_temp ))
		}
	}
	max_acc = max (list_parameter[,3])
	good = list_parameter[which(list_parameter[,3] == max_acc),]
	return (good[1,])
}



SVMLOO = function (data){
	
	
	data_predic = data[,- which(colnames(data)=="drugg")]
	drugg = changeList(data[,which(colnames(data)=="drugg")])
	#data = cbind(data_predic, drugg)
	parameter = optimization (data)

	v_predict = NULL
	for (i_line in seq(dim(data)[1])){
		data_train = data[-i_line,]
		data_test = data_predic[i_line,]
		model_svm = ksvm(drugg~., data_train, kernel="rbfdot", kpar=list(sigma=parameter[1]),C=parameter[2], type = "C-svc")
		y_predicted=predict(model_svm, data_test)
		v_predict = append (v_predict,y_predicted)
	}
	return (changeList(v_predict))
}




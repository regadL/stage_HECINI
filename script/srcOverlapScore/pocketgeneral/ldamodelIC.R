#!/usr/bin/env Rscript


# By BORREL Alexandre
# 05-2012

source ("tool.R")
require(plotrix)
library (rpart)
library (MASS)
source ("elimcor_avecY.R")
source ("predictive_power.R")





cartCumulData = function (res.lda, data_global_train, data_global_test, descriptor){
	
	print (descriptor)	
		
	Y_real_test = changeList (data_global_test$drugg)
	Y_real_train = changeList (data_global_train$drugg)

	Y_real_fusion = append (Y_real_train, Y_real_test)

	data_fusion = rbind (data_global_train, data_global_test)
	data_test_lda = data_fusion [,descriptor]
	data_train_lda = data_global_train [,descriptor]


	data_test_lda = data_test_lda[,-which(colnames(data_test_lda) == "drugg")]
	#data_global_train = data_global_train[,-which(colnames(data_global_train) == "drugg")]
	#data_global_test = data_global_test[,-which(colnames(data_global_test) == "drugg")]

	bad_predicted_train = rownames (data_global_train) [which (res.loo$class !=  Y_real_train)]
	
	v_predict_lda = predict (res.lda, data_test_lda)
	# FUSION
	bad_predicted_test = rownames (data_fusion) [which (v_predict_lda$class !=  Y_real_fusion)]


	data_cart = data_fusion[bad_predicted_test,]
	print (dim (data_cart))

	res.cart = rpart (drugg ~., data = data_cart, control=rpart.control(minsplit=3), method = "class")

	png (paste(data_global_train, "_cart_final.png", sep = ""), 2000, 1700)
	plot (res.cart, main = "Cart model", cex.main = 2, cex = 6)
	text (res.cart, use.n=TRUE, cex = 2.75)
	dev.off()

	#print (which((data_fusion[,"drugg"] == 1)==TRUE))	

	data_drug = data_fusion [which((data_fusion[,"drugg"] == 1)==TRUE),]
	data_drug = na.omit (data_drug)
	data_no_drug = data_fusion [which((data_fusion[,"drugg"] == 0)==TRUE) ,]
	data_no_drug = na.omit (data_no_drug)

	data_drug_bad = data_cart[which((data_cart[,"drugg"] == 1)==TRUE),]
	data_drug_bad = na.omit (data_drug_bad)
	bad_drug_predictic = which((data_cart[,"drugg"] == 1)==TRUE)
	data_drug_good = data_drug[-bad_drug_predictic,]



	data_no_drug_bad = data_cart[which((data_cart[,"drugg"] == 0)==TRUE),]
	bad_predic_no_drug = which((data_cart[,"drugg"] == 0)==TRUE)
	data_no_drug_bad = na.omit (data_no_drug_bad)

	data_no_drug_good = data_no_drug[-bad_predic_no_drug,]
	

	# cart druggable
	bad = append (rep (0,dim(data_drug_good)[1]), rep (1,dim(data_drug_bad)[1]))
	print (bad)
	data_cart_drug = rbind (data_drug_good, data_drug_bad)
	data_cart_drug = cbind (data_cart_drug, bad)
	data_cart_drug = data_cart_drug[,-which(colnames(data_cart_drug) == "drugg")]

	res.cart = rpart (bad ~., data = data_cart_drug, method = "class")

	png (paste(data_global_train, "_cart_druggable_bad.png", sep = ""), 2000, 1700)
	plot (res.cart, main = "", cex.main = 5, cex = 6)
	text (res.cart, use.n=TRUE, cex = 4)
	dev.off()





	# cart no druggable
	bad = append (rep (0,dim(data_no_drug_good)[1]), rep (1,dim(data_no_drug_bad)[1]))
	print (bad)
	data_cart_no_drug = rbind (data_no_drug_good, data_no_drug_bad)
	data_cart_no_drug = cbind (data_cart_no_drug, bad)
	data_cart_no_drug = data_cart_no_drug[,-which(colnames(data_cart_no_drug) == "drugg")]

	res.cart = rpart (bad ~., data = data_cart_no_drug, method = "class")

	png (paste(data_global_train, "_cart_no_druggable_bad.png", sep = ""), 2000, 1700)
	plot (res.cart, main = "", cex.main = 5, cex = 5)
	text (res.cart, use.n=TRUE, cex = 4)
	dev.off()


	# plot avec les descripteurs interet
	#png ("~/Dropbox/Mean_alpha_sphere_SA.png")
	#l = list (data_drug[,"Mean_alpha.sphere_SA"], data_no_drug[,"Mean_alpha.sphere_SA"])
	#multhist(l, col = c(5,2),cex.names = 0.8, freq = FALSE, cex.axis = 1)
	#dev.off()
	#png ("~/Dropbox/positive_residues.png")
	#l = list (data_drug[,"positive_residues"], data_no_drug[,"positive_residues"])
	#multhist(l, col = c(5,2), cex.names = 0.8,freq = FALSE, cex.axis = 1)
	#dev.off()

	#png ("~/Dropbox/nitrogen_drug.png")
	#l = list (data_drug_good[,"nitrogen"], data_no_drug_bad[,"nitrogen"])
	#multhist(l, col = c(5,2), cex.names = 0.8,freq = FALSE, cex.axis = 1)
	#dev.off()
	#png ("~/Dropbox/T_drug.png")
	#l = list (data_drug_good[,"T"], data_no_drug_bad[,"T"])
	#multhist(l, col = c(5,2),cex.names = 0.8, freq = FALSE, cex.axis = 1)
	#dev.off()
	#png ("~/Dropbox/aliphatic_residues_drug.png")
	#l = list (data_drug_good[,"aliphatic_residues"], data_no_drug_bad[,"aliphatic_residues"])
	#multhist(l, col = c(5,2),cex.names = 0.8, freq = FALSE, cex.axis = 1)
	#dev.off()

	#png ("~/Dropbox/F_no_drug.png")
	#l = list (data_no_drug_good[,"F"], data_no_drug_bad[,"F"])
	#multhist(l, col = c(5,2), cex.names = 0.8,freq = FALSE, cex.axis = 1)
	#dev.off()
	#png ("~/Dropbox/Y_no_drug.png")
	#l = list (data_no_drug_good[,"Y"], data_no_drug_bad[,"Y"])
	#multhist(l, col = c(5,2), cex.names = 0.8, freq = FALSE, cex.axis = 1)
	#dev.off()

} 



modelLDA = function (desc, data_test, path_file_test){
	
	Y_out = changeList (desc$drugg)
	desc_lda = desc[,-which(colnames(desc) == "drugg")] # remove drugg column

	res.lda = lda(x=desc_lda,grouping=Y_out) # model
	res.loo = lda(x=desc_lda,grouping=Y_out, CV = TRUE) # leave-one-out
	
	
	bad_predictied = which (res.loo$class !=  Y_out)

	cart_lda_predited (bad_predictied, desc_lda, path_file_test)


	# ACP with 2 dataset
	ACPArea (desc, data_test, res.loo$class, res.lda, path_file_test)

	return (list(res.lda, res.loo))
	
}	


cart_lda_predited = function (bad_predicted, data, path_file_test){
	
	predicted = rep (0,dim (data)[1])
	predicted[bad_predicted] = 1
	data = cbind (data, predicted)
	res.cart = rpart (predicted ~., data = data, method = "class", control=rpart.control(minsplit=5))
	png (paste (path_file_test, "_MODELE_Predicted_CART.png", sep = ""), 2000, 1700)
	plot (res.cart, main = "", cex.main = 2)
	text (res.cart, use.n=TRUE, cex = 1.75)
	dev.off()

}




ACPArea = function (data, data_test, bad_predict, lda_model, path_png){

	bad_predict = which (changeList(data$drugg) != bad_predict)
	# ACP all points
	drug_vector = changeList (data$drugg)
	data = data[,-which(colnames (data)== "drugg")]	
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
	png (paste(path_png, "_dataSet.png", sep = ""), 2480, 3508)
	par(mfrow = c(2,1))	
	par(mar=c(8,8,8,8))


	# With druggable and non-druggable
	plot(data_plot[,1],data_plot[,2], type = "n", main = "", xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75)
	text (data_plot[,1],data_plot[,2], col = color_point, labels = drug_vector, cex = 2.9)	

	#arrows(0, 0, cp[,1]*factor, cp[,2]*factor, col = as.integer(c))
	#text(cp[,1]*factor, cp[,2]*factor, colnames(data.scale), col = as.integer(c), cex = 2.5)
	
	abline(h=0,v=0)
 	#legend("topleft", col=c("blue",2,3,4,6), legend = c("Good Predict", "Bad Predict","desc RADII3","desc perso","desc Fpocket"),pch=c(20,20,26,26,26),lty=c(0,0,1,1,1), cex = 3.5)
	
	# with data test
	plot(data_plot[,1],data_plot[,2], type = "n", main = "", xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75)
	
	drug_vector_test = changeList (data_test$drugg)
	data_test = data_test[,-which(colnames(data_test) == "drugg")]
	y_predict = predict (lda_model, data_test)
	bad_predict = which (drug_vector_test != y_predict$class)
	color_point_test = rep (4, dim (data_test)[1])
	no_drug = which (drug_vector_test == "nd")
	no_drug_bad_predict = intersect (as.vector(no_drug), as.vector(bad_predict))
	color_point_test [bad_predict] = 2
	color_point_test[no_drug_bad_predict] = "#CC06D3"

	
	data_scale_test = scale (data_test)
	data_cor_test=cor(data_scale_test)
	data_eigen_test=eigen(data_cor_test)
	lambda = data_eigen_test$values
	cp = data_eigen_test$vectors
	rownames (cp) = colnames (data_scale_test)
	colnames (cp) = colnames (data_scale_test)
	data_plot = as.matrix(data_scale_test)%*%cp
	
	text (data_plot[,1],data_plot[,2], labels = drug_vector_test, cex = 2.9, col = color_point_test)	
	abline(h=0,v=0)
	dev.off ()

}


modelCART = function(data_train_global, data_train, res.lda, path_file_train){
	Y_out = changeList (data_train$drugg )
	bad_prediction = which (res.lda$class !=  Y_out) # result LDA
	bad_prediction_PDB = rownames (data_train)[bad_prediction]	
	data_bad = data_train_global[bad_prediction_PDB,]
	res.cart = rpart (drugg ~., data = data_bad, control=rpart.control(minsplit=3), method = "class")

	png (paste (path_file_train, "_CORRECTOR.png", sep = ""), 2000, 1700)
	plot (res.cart, main = "Cart model", cex.main = 2, cex = 6)
	text (res.cart, use.n=TRUE, cex = 2.75)
	dev.off()

	return (res.cart) 
}




testModel= function (res.lda, res.cart, data_test, descriptor, print_file, path_file, boundary_inf, boundary_sup){

	Y_real = changeList (data_test$drugg)
	
	color_point = as.matrix(rep ("#1D12BB", dim(data_test)[1]))
	
	rownames (color_point) = rownames (data_test)
	
	data_test_lda = data_test[,descriptor]
	list_drugable = data_test_lda[,which(colnames(data_test_lda) == "drugg")]
	data_test_lda = data_test_lda[,-which(colnames(data_test_lda) == "drugg")] # remove drugg column
	#data_test_without_drugg = data_test[,-which(colnames(data_test) == "drugg")] # remove drugg column
	
	# LDA
	v_predict_lda = predict (res.lda, data_test_lda)


	# AREA ERROR
	data_bad = selectDataCorrected(v_predict_lda, data_test, boundary_inf, boundary_sup)

	

	### write for barplot with name protein
	bad_prediction = which (v_predict_lda$class !=  Y_real)
	write_element = cbind (v_predict_lda$posterior, list_drugable)
	write.table (write_element, file = paste (path_file, "_quality_predict", sep = ""), row.names = TRUE, col.names = TRUE)

	data_bad_real = data_test[bad_prediction,]
	# CART on descriptor  selected
	

	#CART
	v_predict_cart = predict (res.cart, newdata = data_bad[,-which (colnames (data_bad) == "drugg")], type="class")
	predict_cart = changeList (v_predict_cart)

	
	color_point[bad_prediction] = "#00D640"

	if (print_file == 1){
		print ("----------------")
		print ("LDA")
		print ("----------------")
		taux_lda = qualityPredict (list(v_predict_lda$class), list(Y_real))

		print ("---ERROR AREA---")
		print (paste(boundary_inf, "->", boundary_sup, sep = ""))
		print ("----------------")

		# list predicting
		Y_real_area = changeList(data_bad$drugg)
		Y_out_LDA_area = v_predict_lda$class
		names (Y_out_LDA_area) = rownames (data_test)
		Y_out_LDA_area = Y_out_LDA_area[rownames (data_bad)]

		# taux area
		taux_lda = qualityPredict (list(Y_out_LDA_area), list(Y_real_area))
		

		#######################################################
		### Cart difference between bad and good predicting ###
		#######################################################
		cart_lda_predited (bad_prediction, data_test, paste (path_file, "_badVSgood_test_", sep = ""))
	}
		

	if (print_file == 1){
		print ("----------------")
		print ("CART -> ERROR AREA")
		print ("----------------")
	
		# taux for cart
		taux_cart = qualityPredict (list(predict_cart), list(Y_real_area))

		print ("---------------")
		print ("global taux")
		print ("---------------")
		cumulTaux(taux_lda, taux_cart )
	
	}

	# revoir les couleur piour les ACP
	#Y_real_bad = Y_real[bad_prediction] # for color
	#v_bad_predict = v_predict_cart[which(predict_cart != Y_real_bad)]	
	#color_point[names (v_bad_predict),] = "#EB0000"
	return (color_point)
}




selectDataCorrected = function (v_predict_lda, data_test, begin, end){

	list_proba = v_predict_lda$posterior
	zone_test = list_proba[which(list_proba[,1] > begin),]
	zone_test = zone_test[which(zone_test[,1] < end),]
	list_PDB = rownames (zone_test)
	return (data_test[list_PDB,])
}





testCART = function (data_test, res.cart, path_file){

	y_real = changeList (data_test$drugg)
	data_test = data_test[,-which(colnames(data_test) == "drugg")] # remove drugg column
	
	# predict	
	ypred = predict(res.cart, newdata = data_test, type="class")
	ypred = changeList (ypred)
	
	taux_cart = qualityPredict (list(ypred), list(y_real))
}



openSpData = function (path_file_descriptor, data){

	list_descriptor = as.vector(read.table (path_file_descriptor)[,1])
	list_descriptor = append (list_descriptor, "drugg")
	return (data[, list_descriptor])	
}


ACP = function (data, color_point, path_png, col.desc){
	
	drug_vector = data[,which(colnames(data) == "drugg")]
	
	drug_vector = changeList (drug_vector)
	data = data[,-which(colnames(data) == "drugg")] 
	data = scale (data)

	data.cor=cor(data)
	data.eigen=eigen(data.cor)

	lambda = data.eigen$values
	var_cap = cumsum(lambda)/sum(lambda)*100
	
	cp = data.eigen$vectors
	rownames (cp) = colnames (data)
	colnames (cp) = colnames (data)
	data_plot = as.matrix(data)%*%cp
	#print (data_plot)
	factor = factorACP (data_plot, cp)
	png (path_png, 2480, 3508)
	par(mar=c(8,8,8,8))
	par(mfrow = c(2,1))
	plot(data_plot[,1],data_plot[,2], main = "Model Druggability", xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75, type = "n")
	text (data_plot[,1],data_plot[,2], col = color_point, labels = drug_vector, cex = 2.2)	
	legend("topleft", col=c("#1D12BB","#00D640","#EB0000"), legend = c("Good predict LDA", "Good predict CART","Bad predict"),pch=c(19,19,19),lty=c(0,0,0), cex = 3.5)  
	abline(h=0,v=0)	

	plot(data_plot[,1],data_plot[,2], col = color_point, main = "Model Druggability", xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75, type = "n")
	c = as.vector (col.desc[,rownames(cp)])
	arrows(0,0,cp[,1]*factor,cp[,2]*factor, col= as.integer(c))
	text(cp[,1]*factor,cp[,2]*factor, rownames(cp),  col= as.integer(c), cex = 3.5)
	abline(h=0,v=0)
	warnings ()	
	dev.off()
	return ()
}

ACPGlobal = function (data1, data2, color_point, path_png, col.desc){
	
	data = rbind (data1,data2)
	nb_value_sep = dim (data1)[1]

	drug_vector = data[,which(colnames(data) == "drugg")]
	drug_vector = changeList (drug_vector)
	data = data[,-which(colnames(data) == "drugg")] 
	data = scale (data)

	data.cor=cor(data)
	data.eigen=eigen(data.cor)

	lambda = data.eigen$values
	var_cap = cumsum(lambda)/sum(lambda)*100
	
	cp = data.eigen$vectors
	rownames (cp) = colnames (data)
	colnames (cp) = colnames (data)
	data_plot = as.matrix(data)%*%cp
	#print (data_plot)
	factor = factorACP (data_plot, cp)
	png (path_png, 2480, 3508)
	par(mar=c(8,8,8,8))
	par(mfrow = c(2,1))
	plot(data_plot[,1],data_plot[,2], main = "Model Druggability", xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75, type = "n")
	points(data_plot[,1][1:nb_value_sep] + 0.1,data_plot[,2][1:nb_value_sep], pch = 15, cex = 2)
	text (data_plot[,1],data_plot[,2], col = color_point, labels = drug_vector, cex = 2.2)	
	legend("topleft", col=c("#1D12BB","#00D640","#EB0000"), legend = c("Good predict LDA", "Good predict CART","Bad predict"),pch=c(19,19,19),lty=c(0,0,0), cex = 3.5)  
	abline(h=0,v=0)	

	plot(data_plot[,1],data_plot[,2], col = color_point, main = "Model Druggability", xlab = paste("CP1: var=", round(var_cap[1],2),"%",paste=" "), ylab = paste("CP2: var-cum=", round(var_cap[2],2),"%",sep=" "), cex.lab = 4, cex.main = 2, cex.axis = 1.75, type = "n")
	c = as.vector (col.desc[,rownames(cp)])
	arrows(0,0,cp[,1]*factor,cp[,2]*factor, col= as.integer(c))
	text(cp[,1]*factor,cp[,2]*factor, rownames(cp),  col= as.integer(c), cex = 3.5)
	abline(h=0,v=0)
	warnings ()	
	dev.off()

	return ()
}

############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file_train = args[1]
path_file_test = args[2]
path_file_global_train = args[3]
path_file_global_test = args[4]
path_file_test_corected = args[5]
path_file_descriptor = args[6]
inf = args[7]
sup = args[8]
elimcor = args[9]


print (elimcor)


d.col = colorACPByTypeOfDescriptors ()

# open for LDA
data_open = openData (path_file_train, elimcor)
descriptor = data_open[[2]]
data_train = data_open[[1]]
data_test = openData (path_file_test, 0)[[1]]

descriptor = append (descriptor, "drugg")
data_test = data_test[,descriptor]

# open for cart
data_global_train = openData (path_file_global_train, 0)[[1]]
data_global_test = openData (path_file_global_test, 0)[[1]]
data_corrected = openData (path_file_test_corected, 0)[[1]]

# CART with data test to see desccriptor implication between good and bad prediction

res.model = modelLDA(data_train, data_test, path_file_test)
res.model2 = modelLDA(data_train, data_corrected[,colnames (data_train)], paste(path_file_test, "_selected", sep = ""))
res.loo = res.model[[2]]
res.lda = res.model[[1]]

res.cart = modelCART (data_global_train, data_train, res.loo, path_file_train)


#cartCumulData ( res.lda, data_global_train, data_global_test, descriptor) 

color_point_test = testModel (res.lda, res.cart, data_global_test, descriptor, 1, path_file_test, inf, sup)
color_point_train = testModel (res.lda, res.cart, data_global_train, descriptor, 0, path_file_train, inf, sup)

#color_point_corected = testModel (res.lda, res.cart, data_corrected, descriptor, 0, paste (path_file_test, "_corrected", sep = ""), inf, sup)


#ACP global one prediction with all datasets
#ACP (data_global_test, color_point_test, paste (path_file_test, "_test_ACP.png", sep = ""), d.col)
#ACP (rbind(data_global_test, data_global_train), append(color_point_test, color_point_train), paste (path_file_test, "_global_ACP.png", sep = ""), d.col)

#ACP one significative descriptor for bad prediction
#ACP(openSpData (path_file_descriptor, data_global_test), color_point_test, paste (path_file_test, "_hydrophobie.png", sep = ""), d.col)
#ACP(openSpData (path_file_descriptor, data_global_train), color_point_train, paste (path_file_train, "_hydrophobie.png", sep = ""), d.col)
#ACPGlobal(openSpData (path_file_descriptor, data_global_train),openSpData (path_file_descriptor, data_global_test) , append (color_point_train, color_point_test), paste (path_file_train, "_Global_hydrophobie.png", sep = ""), d.col)


#ACPGlobal(openSpData (path_file_descriptor, data_global_test),openSpData (path_file_descriptor, data_corrected) , append (color_point_test,color_point_corected ), paste (path_file_test, "_predicted.png", sep = ""), d.col)

print ("TEST CART ONLY")
print ("--------------")
print ("   --------   ")

#print ("----------CART----------")
#print ("-----TRAINNING DATA-----")
#testCART (data_global_train, res.cart, path_file_train)
#print ("-------TEST DATA--------")
#testCART (data_global_test, res.cart, path_file_test)

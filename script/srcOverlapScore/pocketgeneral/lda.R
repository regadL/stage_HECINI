#!/usr/bin/env Rscript


# By BORREL Alexandre
# 04-2012

library (MASS)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")
source ("learningMachine.R")
require(plotrix)
library ("klaR")



densityData = function (data1, data2, path_png){

	nb_descriptor = dim (data1)[2]
	nb_page = as.integer(nb_descriptor / 9) + 1

	for (page in seq (0, nb_page-1,1)){
		png (paste ( path_png, "_density_", as.character (page), ".png", sep = ""), 2480, 3508)
		par (mfrow = c(4,3))
		for (descriptor in seq (1,12)){
			num_col = page*12 + descriptor
			if (num_col <= nb_descriptor){
				par(mar=c(6,6,6,6))
				xl = calculLimListData (data1[,num_col], data2[,num_col])
				xd = calculLimListData (density(data1[,num_col])$y, density(data2[,num_col])$y)
				plot (density (data1[,num_col]), main = "", xlim = c(xl[1], xl[2]), ylim = c(0, xd[2]), col = 5,lwd=3 , xlab = "", ylab = "", yaxt = "n", xaxt = "n")
				par (new = TRUE)
				plot (density (data2[,num_col]), main = "", xlim = c(xl[1], xl[2]), ylim = c(0, xd[2]), col = 2,lwd=3 , xlab = "", ylab = "", cex.axis = 2)
				title(main=colnames (data1)[num_col], cex.main = 2.5, ylab = "Density fonction", cex.lab = 2)
				legend("topright", col=c(5,2), legend = c("Predict", "Bad Predict"), pch=c(16,16), lty=c(0,0,1,1,1), cex = 2)
			}
		}
		dev.off()
	}


}


barplotDescriptor = function(res.lda, name_file, name_barplot, data_train){

	aa = dim (data_train)
	#print (attributes(res.lda))
	#print ("ggggg")
	#print (res.lda$scaling)
	#print ("dddddd")
	# coeficent LDA
	coef = res.lda$scaling[,1]
				
	# Standardized
	coef = normalizationScalingLDA (coef, data_train)

	vcol = rep(2,length=length(coef))
	names(vcol) = names(coef)
	vcol[which(coef<0)]=4
	coef.sort = sort (abs(coef),decreasing = TRUE)
	print (coef.sort)
	png (paste (name_file, ".png", sep = ""), 2000, 2000)
	par( mar=c(45,10,5,5))
	barplot(coef.sort, names.arg = names(coef.sort), las=2, cex.names = 2.75, col = vcol[names(coef.sort)], main = "", cex.main = 3, cex.axis = 3, ylab = "Significativité des descripteurs (%)", cex.lab = 3)
	#legend("topright",legend=c("coef >0","coef <0"), col=c(2,4),lty=2, cex = 2.75)
	dev.off ()

}


ldaLeaveOneOut = function (desc, data_global, name_png, graph){
	print ("----Leave one out----")
	Y =  as.factor(desc$drugg)
	Y2 = rep("d",length(Y))
	Y2[which(Y==0)]="nd"
	Y_out = Y2

	desc = desc[,-which(colnames(desc) == "drugg")] # remove drugg column
	res.lda = lda(x=desc,grouping=Y2, CV = TRUE) # leave-one-out
	
	print (length (colnames (desc)))
	print (colnames (desc))
	
	# histogram and density
	bad_prediction = which (res.lda$class !=  Y_out)

	

	data_global = data_global[rownames (desc),] # same rown that data predicting	
	
	if (graph == 1){
		histData (data_global[-bad_prediction,], data_global[bad_prediction,], paste(name_png,".pdf", sep = ""))
		#densityData (data_global[-bad_prediction,], data_global[bad_prediction,], name_png)
	
		# ACP global
		ACPPredict (desc, bad_prediction, Y_out, name_png, "ACP LOO Prediction")
	}

	write.table (rownames(data_global[-bad_prediction,]), file = "good_predict", row.names = FALSE, col.names = FALSE)
	write.table (rownames(data_global[bad_prediction,]), file = "bad_predict", row.names = FALSE, col.names = FALSE)

	drug = which(Y_out == "d")
	no_drug =  which(Y_out == "nd")	
	
	write.table (rownames(data_global)[intersect(drug,bad_prediction)], file = "bad_predict_drug", row.names = FALSE, col.names = FALSE)
	write.table (rownames(data_global)[intersect(no_drug,bad_prediction)], file = "bad_predict_no_drug", row.names = FALSE, col.names = FALSE)

	data_good_predict = data_global[-bad_prediction,]
	list_drug = Y_out[-bad_prediction]
	drug_good_predict = which(list_drug == "d")
	no_drug_good_predict = which(list_drug == "nd")

	write.table (rownames(data_good_predict)[drug_good_predict], file = "good_predict_drug", row.names = FALSE, col.names = FALSE)
	write.table (rownames(data_good_predict)[no_drug_good_predict], file = "good_predict_no_drug", row.names = FALSE, col.names = FALSE)

	write_element = cbind (res.lda$posterior, Y)
	write.table (write_element, file = paste (name_png, "_LOO_quality_predict", sep = ""), row.names = TRUE, col.names = TRUE)

	return (list(res.lda$class, Y_out))
}



scrambling = function (data_train, data_test){
	
	Y = as.factor(data_train$drugg)
	Y2 = rep("d",length(Y))
	Y2[which(Y==0)]="nd"
	
	Y1 = as.factor(data_test$drugg)
	Y_out = rep("d",length(Y1))
	Y_out[which(Y1==0)]="nd"
	Y_out = sample (Y_out)	
	
	data_train = data_train[,-which(colnames(data_train) == "drugg")] # remove drugg column
	data_test = data_test[,-which(colnames(data_test) == "drugg")] # remove drugg column
	
	res.lda = lda(x=data_train,grouping=Y2) # LDA model
	v_predict = predict (res.lda, data_test) # prediction

	bad_prediction = which (v_predict$class !=  Y_out)



	return (list(v_predict$class, Y_out, res.lda))
}



ldaTrainTest = function (data_train, data_test, path_filout, name_barplot, draw_plot, name_ACP, graph){
	
	Y = as.factor(data_train$drugg)
	Y2 = rep("d",length(Y))
	Y2[which(Y==0)]="nd"
	
	Y1 = as.factor(data_test$drugg)
	Y_out = rep("d",length(Y1))
	Y_out[which(Y1==0)]="nd"	
	
	data_train_in = data_train[,-which(colnames(data_train) == "drugg")] # remove drugg column
	data_test = data_test[,-which(colnames(data_test) == "drugg")] # remove drugg column
	
	res.lda = lda(x=data_train_in,grouping=Y2) # LDA model
	v_predict = predict (res.lda, data_test) # prediction

	bad_prediction = which (v_predict$class !=  Y_out)

	# if not different between real and predict, graph not possible
	if ( length(bad_prediction) < 5 ){
		graph = 0
	}else{}

	write_element = cbind (v_predict$posterior, Y1)
	write.table (write_element, file = paste (path_filout, "_", name_ACP, "_quality_predict", sep = ""), row.names = TRUE, col.names = TRUE)

	if (graph == 1){
		ACPPredict (data_test, bad_prediction, Y_out, path_filout, paste("Prediction", name_ACP, sep = ""))
		
		# barplot descriptors scaling
 		if (draw_plot == 1){
			barplotDescriptor(res.lda, path_filout, name_barplot, data_train)
		}else{}	
	}else{}

	return (list(v_predict$class, Y_out, res.lda))
}

MDSColor = function (data, list_descriptor, name_file){

	data = data[,-which (colnames(data) == "drugg")]
	MC = cor(data)
	dist1 = 1-MC
	#print (dist1)
	fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
	png (name_file, 2480, 1700)
	par( mar=c(12,12,12,12))
	fit <- cmdscale(as.dist(dist1), eig=TRUE, k=2)
	#print (fit)
	plot (fit$points[,1], fit$points[,2], main="MDS variables selected", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 3, cex.axis = 3, cex.main = 2.75, type = "n")
	text (fit$points[,1], fit$points[,2], rownames (fit$points), cex = 2.8)
	text (fit$points[list_descriptor,1], fit$points[list_descriptor,2], labels = rownames (fit$points)[list_descriptor], cex = 2.8, col = "red")
	dev.off()
}




############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file1 = args[1]
option_leave_one_out = as.integer(args[2])
path_file2 = args[3]
path_filout = args[4]
name_plot = args[5]
elimcor = args[6]
graph = as.integer(args[7])
path_global = args[8]

data_open = openData (path_file1, elimcor)

desc = data_open[[1]]
descriptor = data_open[[2]]

data_global = read.table (path_global, header=TRUE, sep = "\t")

if (option_leave_one_out == 1){
	print (paste ("elimCor -> ", elimcor, sep = ""))
	out_lda = ldaLeaveOneOut (desc, data_global, path_filout, graph)
	out_del = qualityPredict (out_lda[1], out_lda[2])

}else{
	desc_test = openData (path_file2, 0)[[1]]
	if (elimcor != 0){
		if (graph == 1){
			MDSColor (openData (path_file1, 0)[[1]], descriptor, paste (path_file1, "_MDS.png"))
		}		
		drugg = desc_test$drugg
		desc_test = desc_test[,descriptor]
		desc_test = cbind (desc_test, drugg)
	}
	

	print ("--------------------------------")
	print ("----Data test and data train----")
	print ("******TRAIN DATASET*****")
	out_lda = ldaTrainTest (desc, desc, paste(path_filout, "train", sep = ""), name_plot, 1, "Train", graph)
	qualityPredict (out_lda[1], out_lda[2])

	print ("--------------------------------")
	print ("*****SAMPLING TRAIN DATASET*****")
	out_trainlda_scr = scrambling (desc, desc)
	out_traindel_scr = qualityPredict (out_trainlda_scr[1], out_trainlda_scr[2])	



	print ("**********TEST DATASET*********")
	out_lda = ldaTrainTest (desc, desc_test, paste(path_filout, "test", sep = ""), name_plot, 0, "Test", graph)
	out_del = qualityPredict (out_lda[1], out_lda[2])

	print ("--------------------------------")
	print ("*****SAMPLING TEST DATASET******")
	out_lda_scr = scrambling (desc, desc_test)
	out_del_scr = qualityPredict (out_lda_scr[1], out_lda_scr[2])


	# save model LDA for use 
	model.lda = out_lda[3]
	
	# plot separation data with descriptors

	desc$drugg[which(desc$drugg == 0)] = "nd"
	desc$drugg[which(desc$drugg == 1)] = "d"

	#print (desc$drugg == 0)
	desc$drugg = as.factor (desc$drugg)
	#print (desc$drugg)
	png (paste(path_file1, "_separation_desc", ".png", sep = ""), 1700, 1500)
	partimat(drugg ~ ., data = desc, method = "lda", plot.matrix = TRUE, imageplot = FALSE, cex = 2, cex.lab = 10, cex.name = 50)
	#drawparti(drugg~., method = "lda", prec = 100, xlab = NULL, ylab = NULL, col.correct = "black", col.wrong = "red", col.mean = "black", col.contour = "darkgrey", gs = as.character(grouping), pch.mean = 19, cex.mean = 1.3, print.err = 0.7, legend.err = FALSE, legend.bg = "white",  imageplot = TRUE, image.colors = cm.colors(nc)) 
	#drawparti (desc$drugg, desc[,2], desc[,3], method = "lda", xlab = colnames(desc)[2], ylab = colnames(desc)[3], pch.mean = 19, cex = 4 )#,imageplot== TRUE)#, image.colors = "white")

	dev.off ()
	save(model.lda, file = paste(path_filout, ".Rdata", sep = ""))
}
#warnings()



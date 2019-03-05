#!/usr/bin/env Rscript


# By BORREL Alexandre
# 07-2012

library (MASS)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")



generateEveryModel4 = function(data_train, data_test, data_global){

	# real LOO
	class = as.factor(data_global$drugg)
	class_real_global = rep("d",length(class))
	class_real_global[which(class==0)]="nd"

	# real train	
	class = as.factor(data_train$drugg)
	class_real_train = rep("d",length(class))
	class_real_train[which(class==0)]="nd"

	# real test
	class = as.factor(data_test$drugg)
	class_real_test = rep("d",length(class))
	class_real_test[which(class==0)]="nd"

	
	list_descriptor = colnames (data_train)[-which(colnames (data_train) == "drugg")] # list descriptor without class descriptor
	nb_descriptor = length (list_descriptor) - 1	


	i_des1 = 1
	while (i_des1 <= nb_descriptor){
		i_des2 = i_des1 + 1
		while (i_des2 <= nb_descriptor){
			i_des3 = i_des2 + 1
			while (i_des3 <= nb_descriptor){
				i_des4 = i_des3 + 1
				while (i_des4 <= nb_descriptor){

					list_des_model = c(list_descriptor[i_des1],list_descriptor[i_des2],list_descriptor[i_des3], list_descriptor[i_des4])
					
					# model
					model.lda = lda(x=data_train[, list_des_model],grouping=class_real_train) # LDA model
					v_predict_train = predict (model.lda, data_train[, list_des_model]) # prediction on training set
					v_predict_test = predict (model.lda, data_test[, list_des_model]) # prediction on training set
					
					# coeficent LDA
					coef = model.lda$scaling[,1]
				
					# Standardized
					coef = normalizationScalingLDA (coef, data_train)


					# rate train
					rate_predict_train = calculTaux2 (v_predict_train$class, class_real_train)	
					acc_train = accuracy(rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])
					se_train = sensibility(rate_predict_train[1], rate_predict_train[4])
					sp_train = sensibility(rate_predict_train[2], rate_predict_train[3])
					MCC_train = MCC  (rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])

					# rate loo
					model.lda.loo = lda(x=data_global[, list_des_model],grouping=class_real_global, CV = TRUE) # leave-one-out
		
					rate_loo = calculTaux2 (model.lda.loo$class, class_real_train)
					acc_loo = accuracy(rate_loo[1], rate_loo[2], rate_loo[3], rate_loo[4])
					se_loo = sensibility(rate_loo[1], rate_loo[4])
					sp_loo = sensibility(rate_loo[2], rate_loo[3])
					MCC_loo = MCC(rate_loo[1], rate_loo[2], rate_loo[3], rate_loo[4])
		
					# rate test
					rate_predict_test = calculTaux2 (v_predict_test$class, class_real_test)	
					acc_test = accuracy(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])	
					se_test = sensibility(rate_predict_test[1], rate_predict_test[4])
					sp_test = sensibility(rate_predict_test[2], rate_predict_test[3])
					MCC_test = MCC(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])	

					# retrieve condition
					#if (acc_loo > 0.82 && acc_train > 0.82 && acc_test > 0.82){
						print ("descriptor")
						print (list_des_model)
						print (as.vector(abs(coef[list_des_model])))
						print ("Acc_loo --- Acc_train --- Acc_test")
						print (paste (acc_loo, acc_train, acc_test, sep = "---"))
						print ("Se_loo --- Sp_loo --- Se_train --- Sp_train --- Se_test --- Sp_test")
						print (paste (se_loo, sp_loo, se_train, sp_train, se_test, sp_test, sep = "---"))
						print ("MCC_loo --- MCC_train --- MCC_test")
						print (paste (MCC_loo, MCC_train, MCC_test, sep = "---"))
						print ("**********************************************************************")
					#}


					i_des4 = i_des4 + 1
				}
				i_des3 = i_des3 + 1
			}
			i_des2 = i_des2 + 1
		}
		i_des1 = i_des1 + 1
	}


}


generateEveryModel5 = function(data_train, data_test, data_global){

	# real LOO
	class = as.factor(data_global$drugg)
	class_real_global = rep("d",length(class))
	class_real_global[which(class==0)]="nd"

	# real train	
	class = as.factor(data_train$drugg)
	class_real_train = rep("d",length(class))
	class_real_train[which(class==0)]="nd"

	# real test
	class = as.factor(data_test$drugg)
	class_real_test = rep("d",length(class))
	class_real_test[which(class==0)]="nd"

	
	list_descriptor = colnames (data_train)[-which(colnames (data_train) == "drugg")] # list descriptor without class descriptor
	nb_descriptor = length (list_descriptor) - 1	


	i_des1 = 1
	while (i_des1 <= nb_descriptor){
		i_des2 = i_des1 + 1
		while (i_des2 <= nb_descriptor){
			i_des3 = i_des2 + 1
			while (i_des3 <= nb_descriptor){
				i_des4 = i_des3 + 1
				while (i_des4 <= nb_descriptor){
					i_des5 = i_des4 + 1
					while (i_des5 <= nb_descriptor){

						list_des_model = c(list_descriptor[i_des1],list_descriptor[i_des2],list_descriptor[i_des3], list_descriptor[i_des4], list_descriptor[i_des5])
						
						# model
						model.lda = lda(x=data_train[, list_des_model],grouping=class_real_train) # LDA model
						v_predict_train = predict (model.lda, data_train[, list_des_model]) # prediction on training set
						v_predict_test = predict (model.lda, data_test[, list_des_model]) # prediction on training set
					
						# coeficent LDA
						coef = model.lda$scaling[,1]
				
						# Standardized
						coef = normalizationScalingLDA (coef, data_train)


						# rate train
						rate_predict_train = calculTaux2 (v_predict_train$class, class_real_train)	
						acc_train = accuracy(rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])
						se_train = sensibility(rate_predict_train[1], rate_predict_train[4])
						sp_train = sensibility(rate_predict_train[2], rate_predict_train[3])
						MCC_train = MCC  (rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])

						# rate loo
						model.lda.loo = lda(x=data_global[, list_des_model],grouping=class_real_global, CV = TRUE) # leave-one-out
		
						rate_loo = calculTaux2 (model.lda.loo$class, class_real_train)
						acc_loo = accuracy(rate_loo[1], rate_loo[2], rate_loo[3], rate_loo[4])
						se_loo = sensibility(rate_loo[1], rate_loo[4])
						sp_loo = sensibility(rate_loo[2], rate_loo[3])
						MCC_loo = MCC(rate_loo[1], rate_loo[2], rate_loo[3], rate_loo[4])
		
						# rate test
						rate_predict_test = calculTaux2 (v_predict_test$class, class_real_test)	
						acc_test = accuracy(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])	
						se_test = sensibility(rate_predict_test[1], rate_predict_test[4])
						sp_test = sensibility(rate_predict_test[2], rate_predict_test[3])
						MCC_test = MCC(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])	

						# retrieve condition
						#if (acc_loo > 0.82 && acc_train > 0.82 && acc_test > 0.82){
						print ("descriptor")
						print (list_des_model)
						print (as.vector(abs(coef[list_des_model])))
						print ("Acc_loo --- Acc_train --- Acc_test")
						print (paste (acc_loo, acc_train, acc_test, sep = "---"))
						print ("Se_loo --- Sp_loo --- Se_train --- Sp_train --- Se_test --- Sp_test")
						print (paste (se_loo, sp_loo, se_train, sp_train, se_test, sp_test, sep = "---"))
						print ("MCC_loo --- MCC_train --- MCC_test")
						print (paste (MCC_loo, MCC_train, MCC_test, sep = "---"))
						print ("**********************************************************************")
						#}

						i_des5 = i_des5 + 1
					}
					i_des4 = i_des4 + 1
				}
				i_des3 = i_des3 + 1
			}
			i_des2 = i_des2 + 1
		}
		i_des1 = i_des1 + 1
	}


}



generateEveryModel3 = function(data_train, data_test, data_global){

	# real LOO
	class = as.factor(data_global$drugg)
	class_real_global = rep("d",length(class))
	class_real_global[which(class==0)]="nd"

	# real train	
	class = as.factor(data_train$drugg)
	class_real_train = rep("d",length(class))
	class_real_train[which(class==0)]="nd"

	# real test
	class = as.factor(data_test$drugg)
	class_real_test = rep("d",length(class))
	class_real_test[which(class==0)]="nd"

	
	list_descriptor = colnames (data_train)[-which(colnames (data_train) == "drugg")] # list descriptor without class descriptor
	nb_descriptor = length (list_descriptor) - 1	


	i_des1 = 1
	while (i_des1 <= nb_descriptor){
		i_des2 = i_des1 + 1
		while (i_des2 <= nb_descriptor){
			i_des3 = i_des2 + 1
			while (i_des3 <= nb_descriptor){

				list_des_model = c(list_descriptor[i_des1],list_descriptor[i_des2],list_descriptor[i_des3])
					
				# model
				model.lda = lda(x=data_train[, list_des_model],grouping=class_real_train) # LDA model
				v_predict_train = predict (model.lda, data_train[, list_des_model]) # prediction on training set
				v_predict_test = predict (model.lda, data_test[, list_des_model]) # prediction on training set
				
				# coeficent LDA
				coef = model.lda$scaling[,1]
				
				# Standardized
				coef = normalizationScalingLDA (coef, data_train)
				
				# Standardized
				#sd_value = apply (data_train[,1:(dim(data_train)[2])], 2, sd)
				#print (sd_value)
				#coef = coef * sd_value
	
				# rate train
				rate_predict_train = calculTaux2 (v_predict_train$class, class_real_train)	
				acc_train = accuracy(rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])
				se_train = sensibility(rate_predict_train[1], rate_predict_train[4])
				sp_train = sensibility(rate_predict_train[2], rate_predict_train[3])
				MCC_train = MCC  (rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])

				# rate loo
				model.lda.loo = lda(x=data_global[, list_des_model],grouping=class_real_global, CV = TRUE) # leave-one-out
		
				rate_loo = calculTaux2 (model.lda.loo$class, class_real_train)
				acc_loo = accuracy(rate_loo[1], rate_loo[2], rate_loo[3], rate_loo[4])
				se_loo = sensibility(rate_loo[1], rate_loo[4])
				sp_loo = sensibility(rate_loo[2], rate_loo[3])
				MCC_loo = MCC(rate_loo[1], rate_loo[2], rate_loo[3], rate_loo[4])
		
				# rate test
				rate_predict_test = calculTaux2 (v_predict_test$class, class_real_test)	
				acc_test = accuracy(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])	
				se_test = sensibility(rate_predict_test[1], rate_predict_test[4])
				sp_test = sensibility(rate_predict_test[2], rate_predict_test[3])
				MCC_test = MCC(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])	

				# retrieve condition
				#if (acc_loo > 0.82 && acc_train > 0.82 && acc_test > 0.82){
					print ("descriptor")
					print (list_des_model)
					print (as.vector(abs(coef[list_des_model])))
					print ("Acc_loo --- Acc_train --- Acc_test")
					print (paste (acc_loo, acc_train, acc_test, sep = "---"))
					print ("Se_loo --- Sp_loo --- Se_train --- Sp_train --- Se_test --- Sp_test")
					print (paste (se_loo, sp_loo, se_train, sp_train, se_test, sp_test, sep = "---"))
					print ("MCC_loo --- MCC_train --- MCC_test")
					print (paste (MCC_loo, MCC_train, MCC_test, sep = "---"))
					print ("**********************************************************************")
				#}
				i_des3 = i_des3 + 1
			}
			i_des2 = i_des2 + 1
		}
		i_des1 = i_des1 + 1
	}
}







modelLDAStepWise = function (data_train, data_test, list_descriptor_selected){ # generated one best model
	# model one one append descriptor
	list_value_sum = NULL
	list_acc_loo = NULL
	list_acc_train = NULL
	list_acc_test = NULL
	list_sp_loo = NULL
	list_sp_train = NULL
	list_sp_test = NULL
	list_se_loo = NULL
	list_se_train = NULL
	list_se_test = NULL

	
	# real train	
	class = as.factor(data_train$drugg)
	class_real_train = rep("d",length(class))
	class_real_train[which(class==0)]="nd"

	# real test
	class = as.factor(data_test$drugg)
	class_real_test = rep("d",length(class))
	class_real_test[which(class==0)]="nd"	


	data_train = data_train[,-which(colnames(data_train) == "drugg")] # remove drugg column for train dataset
	data_test = data_test[,-which(colnames(data_test) == "drugg")] # remove drugg column for train dataset
	list_descriptors = colnames (data_train) # without class and without first descriptor (need 2 descriptors for LDA model)
	
	for (desc_selected in list_descriptor_selected){
		list_descriptors = list_descriptors[-which(list_descriptors == desc_selected)]
		
	}

	for (descriptor in list_descriptors){
		list_descriptor = append(list_descriptor_selected, descriptor)
		
		data_temp_test = data_test[,list_descriptor]
		data_temp_train = data_train[,list_descriptor]
		res.lda = lda(x=data_temp_train,grouping=class_real_train) # LDA model
		v_predict_train = predict (res.lda, data_temp_train) # prediction on training set
		v_predict_test = predict (res.lda, data_temp_test) # prediction on training set

		# rate train
		rate_predict_train = calculTaux2 (v_predict_train$class, class_real_train)	
		acc_train = accuracy(rate_predict_train[1], rate_predict_train[2], rate_predict_train[3], rate_predict_train[4])
		se_train = sensibility(rate_predict_train[1], rate_predict_train[4])
		sp_train = sensibility(rate_predict_train[2], rate_predict_train[3])

		# rate loo
		res.lda.loo = lda(x=data_temp_train,grouping=class_real_train, CV = TRUE) # leave-one-out
		
		rate_loo = calculTaux2 (res.lda.loo$class, class_real_train)
		acc_loo = accuracy(rate_loo[1], rate_loo[2], rate_loo[3], rate_loo[4])
		se_loo = sensibility(rate_loo[1], rate_loo[4])
		sp_loo = sensibility(rate_loo[2], rate_loo[3])
		
		# rate test
		rate_predict_test = calculTaux2 (v_predict_test$class, class_real_test)	
		acc_test = accuracy(rate_predict_test[1], rate_predict_test[2], rate_predict_test[3], rate_predict_test[4])	
		se_test = sensibility(rate_predict_test[1], rate_predict_test[4])
		sp_test = sensibility(rate_predict_test[2], rate_predict_test[3])

	

		# list to select best append descriptor		
		list_value_sum = append (list_value_sum, acc_loo + acc_train + acc_test*0.6)
		
		#loo
		list_acc_loo =  append (list_acc_loo, acc_loo)
		list_sp_loo = append (list_sp_loo, sp_loo)
		list_se_loo = append (list_se_loo, se_loo)
		#train
		list_acc_train =  append (list_acc_train, acc_train)
		list_sp_train = append (list_sp_train, sp_train)
		list_se_train = append (list_se_train, se_train)
		#test
		list_acc_test = append (list_acc_test, acc_test)
		list_sp_test = append (list_sp_test, sp_test)
		list_se_test = append (list_se_test, se_test)

	}
	list_train_loo = list_acc_loo + 0.8*list_acc_train # change selection criteria
	i_max_acc = which (list_train_loo == max(list_train_loo))[1]
	return (list(list_descriptors[i_max_acc], list_acc_loo[i_max_acc], list_acc_train[i_max_acc], list_acc_test[i_max_acc],list_sp_loo[i_max_acc],list_sp_train[i_max_acc],list_sp_test[i_max_acc],list_se_loo[i_max_acc],list_se_train[i_max_acc],list_se_test[i_max_acc]))
}




StepWise =function (data_train, data_test){

	nb_descriptor = dim(data_train)[2] - 2

	# Append descriptor more significative

	data_train_group1 = separeData (data_train, "drugg")[[1]]
	data_train_group2 = separeData (data_train, "drugg")[[2]]
	desc = moreSignif (data_train_group1, data_train_group2)

	list_descriptor = desc
	acc_LOO = NULL
	acc_TRAIN = NULL
	acc_TEST = NULL
	sp_LOO = NULL
	sp_TRAIN = NULL
	sp_TEST = NULL
	se_LOO = NULL
	se_TRAIN = NULL
	se_TEST = NULL


	while (nb_descriptor >=1 ){
		result_desc_selected = selectBestDescriptorStepWise (data_train, data_test, list_descriptor)
	
		# append value
		list_descriptor = append (list_descriptor, result_desc_selected[[1]])
	
		acc_LOO = append (acc_LOO, result_desc_selected[[2]])
		acc_TRAIN = append (acc_TRAIN, result_desc_selected[[3]])
		acc_TEST =  append (acc_TEST, result_desc_selected[[4]])

		sp_LOO = append (sp_LOO, result_desc_selected[[5]])
		sp_TRAIN = append (sp_TRAIN, result_desc_selected[[6]])
		sp_TEST =  append (sp_TEST, result_desc_selected[[7]])

		se_LOO = append (se_LOO, result_desc_selected[[8]])
		se_TRAIN = append (se_TRAIN, result_desc_selected[[9]])
		se_TEST =  append (se_TEST, result_desc_selected[[10]])

		nb_descriptor = nb_descriptor - 1

	}


	png (paste(path_file_train, "_selectACC_descriptor.png"), 2480, 1700)
	par (xaxs="i",yaxs="i")

	plot (seq(2,length(acc_LOO) +1), acc_LOO, col = "red", ylim = c(0,1), xlim = c(0, length(acc_LOO) + 2), xlab = "Descriptors", ylab = "%", cex.axis = 2.5, cex.lab = 2.5, pch = 19, cex = 2)

	points (seq(2,length(acc_LOO) +1), acc_TRAIN, col = "orange", pch = 19, cex = 2)
	points(seq(2,length(acc_LOO) +1), acc_TEST, col = "blue", pch = 19, cex = 2)
	legend ("topleft", legend = c("Accuracy LOO", "Accuracy TRAIN", "Accuracy TEST"),col = c("red", "orange", "blue"), cex = 2.75, bty="y", pt.cex = 2, pch = 19)
	grid(length(acc_LOO) + 2, 10, col = 1, lwd = 2) 
	graphics.off()


	png (paste(path_file_train, "_selectSPSE_descriptor.png"), 2480, 1700)
	par (xaxs="i",yaxs="i", mfrow = c(2,2))

	plot (seq(2,length(sp_LOO) + 1), sp_LOO, pch = 19, col = "red", ylim = c(0,1), xlim = c(0, length(acc_LOO) + 2), xlab = "NB descriptor", ylab = "%", cex.lab = 2.5, cex.axis = 2.5, lwd = 4, main = "LOO")
	points (seq(2,length(se_LOO) + 1), se_LOO, col = "red", pch = 8, cex = 2)
	grid(length(sp_LOO) + 2, 10, col = 1, lwd = 2) 

	plot (seq(2,length(sp_TRAIN) + 1), sp_LOO, pch = 19, col = "orange", ylim = c(0,1), xlim = c(0, length(acc_LOO) + 2), xlab = "NB descriptor", ylab = "%", cex.lab = 2.5, cex.axis = 2.5, lwd = 4, main = "TRAIN")
	points (seq(2,length(se_TRAIN) + 1), se_TRAIN, col = "orange", pch = 8, cex = 2)
	grid(length(sp_LOO) + 2, 10, col = 1, lwd = 2) 

	plot (seq(2,length(sp_LOO) + 1), sp_TEST, pch = 19, col = "blue", ylim = c(0,1), xlim = c(0, length(acc_LOO) + 2), xlab = "NB descriptor", ylab = "%", cex.lab = 2.5, cex.axis = 2.5, lwd = 4, main = "TEST")
	points (seq(2,length(se_TEST) + 1), se_TEST, col = "blue", pch = 8, cex = 2)
	grid(length(sp_LOO) + 2, 10, col = 1, lwd = 2) 

	plot (100,100,type = "n")
	legend ("topleft", legend = "SP",col = "black", cex = 2.75, bty="y", pt.cex = 2, pch = 19)
	legend ("topright", legend = "SE",col = "black", cex = 2.75, bty="y", pt.cex = 2, pch = 8)
	graphics.off()




	list_sum = acc_LOO + acc_TRAIN # max sum acc on loo and training set
	i_best = which(list_sum == max(list_sum) )[1] # first element for minimize number of descriptors (!!!!!)

	# case only one descriptor
	if (i_best == 1){
		i_best = 2
	}


	if (nb_descriptor_selected != 0){
		list_out = list_descriptor[seq (1,nb_descriptor_selected)]
		print (list_out)
	}else{
		list_out = list_descriptor[seq (1,i_best)]
		print (list_out)
	}

	print ("accuracy")
	print (paste (acc_LOO[i_best], acc_TRAIN[i_best], acc_TEST[i_best], sep = " "))

	# MDS -> visualisation descriptor position on MDS

	MDSDesc(data_global, list_out, descriptor, paste(path_file_train,"_MDS.png"))




}



MDSDesc = function (data_global, list_descriptor_selected, descriptor_in, name_file){

	data_global = data_global[,-which (colnames(data_global) == "drugg")]
	MC = cor(data_global)
	dist1 = 1-MC
	fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
	png (name_file, 2480, 1700)
	par( mar=c(12,12,12,12))
	fit <- cmdscale(as.dist(dist1), eig=TRUE, k=2)
	plot (fit$points[,1], fit$points[,2], main="", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 3, cex.axis = 3, cex.main = 2.75, type = "n")
	text (fit$points[,1], fit$points[,2], rownames (fit$points), cex = 2.8)
	
	#print (des_selected)
	#print (descriptor_in)
	#print (fit$points[which(rownames(fit$points) == des_selected),])


	text (fit$points[descriptor_in,1], fit$points[descriptor_in,2],descriptor_in, cex = 2.8, col = "blue")
	print (descriptor_in)
	print (fit$points[list_descriptor_selected,2])
	text (fit$points[list_descriptor_selected,1], fit$points[list_descriptor_selected,2], list_descriptor_selected, cex = 2.8, col = "red")


	dev.off()
}





###########
## MAIN  ##
###########

args <- commandArgs(TRUE)
path_file_train = args[1]
path_file_test = args[2]
option_selected = args[3]
#nb_descriptor_selected = args[4] a voir peut etre plus tard


##############
# Open Files #
##############

data_train_open = openData (path_file_train, 0)

#print (data_train_open)

data_train = data_train_open[[1]]
descriptor = data_train_open[[2]]
data_test = openData (path_file_test,0)[[1]]

print (descriptor)

data_global = rbind (data_train, data_test)


print ("For control")
print(dim (data_global))


########
# RUN  #
########

if (option_selected == "AllModel"){
	generateEveryModel5 (data_train, data_test, data_global)
}else{
	StepWise (data_train, data_test)
}






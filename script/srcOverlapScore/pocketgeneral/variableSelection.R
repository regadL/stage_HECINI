#!/usr/bin/env Rscript

# By BORREL Alexandre
# 07-2012

source ("tool.R")
source ("predictive_power.R")
source ("learningMachine.R")


# faire en meme temps que la selection un MDS avec les differentes data



# retrieve descriptor more significant
moreSignif = function (data1, data2){
	
	p_val = 100
	des_out = NULL

	nb_descriptor = dim (data2)[2]
	for (i_des in seq (1,nb_descriptor)){
		if (conditionTtest(data1[,i_des], data2[,i_des]) == 1){
			pval_temp = comparisonTest(data1[,i_des], data2[,i_des], "parametric")
		}else if (conditionTtest(data1[,i_des], data2[,i_des]) == 0){
			pval_temp = comparisonTest(data1[,i_des], data2[,i_des], "non-parametric")
		}else{
			pval_temp = 100
		}
		if (pval_temp < p_val){
			des_out = i_des
			p_val = pval_temp
		}
	}
	return (colnames (data1)[des_out])
}

# retrieve descriptor less correlated 
lessCor = function (data, list_descriptor, class_group){
	data = data[,-which(colnames (data) == class_group)]
	cor = 100
	des_out = NULL

	nb_descriptor = dim (data)[2]
	for (i_des in seq (1,nb_descriptor)){
		cor_temp = 0
		for (descriptor in list_descriptor){
			cor_temp = cor_temp + abs(cor.test (data[,which(colnames (data) == descriptor)], data[,i_des])$estimate)
		}
		
		if (cor_temp <= cor){
			cor = cor_temp
			des_out = i_des
		}
	}
	return (colnames (data)[des_out])
}


# retieve descriptor correlated with descriptor
descriptorCorreleted = function (data, descriptor, value_cor, class_des){
	list_des_out = NULL
	nb_descriptor = dim (data)[2]
	
	for (i_des in seq (1,nb_descriptor)){
		if (descriptor == colnames(data)[i_des] || class_des == colnames(data)[i_des]){
			next
		}
		cor_value = cor.test (data[,descriptor], data[,i_des])
		if (abs(cor_value$estimate) >= value_cor){
			list_des_out = append (list_des_out, colnames (data)[i_des])
		}
	}
	list_des_out = append (list_des_out, descriptor)
	return (list_des_out)
}



# search best combination
combination = function (data_train, data_test, list_des_init, list_desc, class_descriptor){

	best_acc = 0
	data_model_train = data_train[,list_des_init]
	data_model_test = data_test[,list_des_init]
	

	for (desc_test in list_desc){
		# train
		data_model_train = data_train[,list_des_init]
		data_model_train = cbind (data_model_train, data_train[,desc_test])
		data_model_train = cbind (data_model_train, data_train[,class_descriptor])
		colnames (data_model_train) = c(list_des_init,desc_test, class_descriptor)

		# test
		data_model_test = data_test[,list_des_init]
		data_model_test = cbind (data_model_test, data_test[,desc_test])
		data_model_test = cbind (data_model_test, data_test[,class_descriptor])
		colnames (data_model_test) = c(list_des_init, desc_test, class_descriptor)

		acc_temp = LDAACC (data_model_train, data_model_test, class_descriptor)
		sum_acc_temp = acc_temp[1] + acc_temp[2]
		if (sum_acc_temp >= best_acc){
			desc_out = desc_test
		}
	}
	return (desc_out)
}




########
# MAIN #
########

args <- commandArgs(TRUE)
path_file_train = args[1]
path_file_test = args[2]
cor_group = as.double(args[3])
class_desc = args[4]


# Open data
data_train = openData(path_file_train, 0)[[1]]
#print (dim(data_train))
data_test = openData(path_file_test, 0)[[1]]
#print (dim(data_test))

# Separe data
data_train_group1 = separeData (data_train, class_desc)[[1]]
data_train_group2 = separeData (data_train, class_desc)[[2]]



#data_train_group1 = data_train_group1[,-which(colnames(data_train_group1) == class_desc)]
#data_train_group2 = data_train_group2[,-which(colnames(data_train_group2) == class_desc)]

# select more significant
desc_init = moreSignif (data_train_group1, data_train_group2)
list_desc_out = NULL
list_desc_out = append (list_desc_out, desc_init)

nb_descriptor = dim (data_train)[2] - 1


while (nb_descriptor > 1 ){

	# Separe data
	data_train_global = data_train
	for (des_remove in list_desc_out){
		data_train_global = data_train_global[,-which(colnames(data_train_global) == des_remove)]
	}
	#print ("*******")
	#print (dim(data_train_global))
	#print ("*******")
	data_train_group1 = separeData (data_train_global, class_desc)[[1]]
	data_train_group2 = separeData (data_train_global, class_desc)[[2]]
	desc = moreSignif (data_train_group1, data_train_group2)

	# remove correled variable
	flag = 0
	for (des_selected in list_desc_out){
		#print (des_selected)
		coreleted = abs(cor.test (data_train[,which(colnames (data_train) == des_selected)], data_train[,which(colnames (data_train) == desc)])$estimate)
		if (coreleted > 0.7){
			flag = 1
		}
	}
	if (flag == 1){
		data_train = delCol (data_train, desc)
		nb_descriptor = nb_descriptor - 1
		next
	}


	#print ("-------")	
	#print (desc)
	#print ("-------")
	# group of descriptor
	list_corr_descriptor = descriptorCorreleted (data_train, desc, cor_group, class_desc)
	#print ("1111")	
	#print (list_corr_descriptor)
	#print ("1111")
	# retrieve best descriptor in group	
	conserved_descriptor = combination (data_train, data_test, list_desc_out, list_corr_descriptor, class_desc)
	#print ("3333333")	
	#print (conserved_descriptor)
	#print ("3333333")	

	list_desc_out = append (list_desc_out, conserved_descriptor)
	if (length (list_corr_descriptor) == 1){
		nb_descriptor = nb_descriptor - length (list_corr_descriptor)
		next
	}else{
		list_remove_descriptor = list_corr_descriptor[-which(conserved_descriptor == list_corr_descriptor)]
		list_intersect = intersect (list_remove_descriptor, list_desc_out)

		# interesection betwwen previous descriptor remove list
		for (inter in list_intersect){

			#print ("888888")			
			#print (inter)
			#print (list_remove_descriptor)
			#print ("888888")
			#print (which(list_remove_descriptor == inter))
			
			list_remove_descriptor = list_remove_descriptor[-which(list_remove_descriptor == inter )]
		}
	}

	nb_descriptor = nb_descriptor - length (list_corr_descriptor)
	
	# del correlated descriptor
	data_train = delCol (data_train, list_remove_descriptor)
	#print ("------")
	#print (list_remove_descriptor)
	#print (dim (data_train))
	#print ("-----")
}

list_desc_out = append (list_desc_out, class_desc)

data_train_LDA = data_train[,list_desc_out]
data_test_LDA = data_test[,list_desc_out]

acc_temp = LDAACC (data_train_LDA, data_test_LDA, class_desc)
acc_loo = LDALOO (data_train_LDA, class_desc)

print (list_desc_out)

print ("accuracy")
print (acc_temp)
print (acc_loo)

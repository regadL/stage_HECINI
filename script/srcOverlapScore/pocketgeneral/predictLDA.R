#!/usr/bin/env Rscript


# By BORREL Alexandre
# 09-2012

library (MASS)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")
source ("learningMachine.R")




############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file_model = args[1]
path_file_descriptor = args[2]
graphics_option = args[3]

print (path_file_model)
print (path_file_descriptor)

# load environment
load(path_file_model,.GlobalEnv)


list_descriptor = try(rownames(model.lda[[1]]$scaling))
list_descriptor = try (rownames(model.lda$scaling))

print (list_descriptor)


# load descriptors
data_apo = openDataDescriptor (path_file_descriptor, list_descriptor)[[1]]

if (dim (data_apo)[1] == 1){
	print ("----proba-----")
	print (predict (model.lda[[1]], data_apo[, list_descriptor])$posterior)	

}else {

	# Y real
	Y =  as.factor(data_apo$drugg)
	Y2 = rep("d",length(Y))
	Y2[which(Y==0)]="nd"
	Y_real = Y2
	
	# remove col drug
	data_apo_predict = data_apo[,-which(colnames(data_apo) == "drugg")] # remove drugg column
	
	#print (data_apo_predict)	

	# prediction
	v_predict = try (predict (model.lda[[1]], newdata = data_apo_predict))
	v_predict = try (predict (model.lda, newdata = data_apo_predict))
	
	# quality
	qualityPredictList (v_predict$class, Y_real)
	
	# bad predicted
	if (graphics_option == 1){ 
		# proba files
		#print (data_apo[rownames (v_predict$posterior), "drugg"] + 1)
		write.table (cbind(v_predict$posterior,data_apo[rownames (v_predict$posterior), "drugg"] + 1), file = paste(path_file_descriptor, ".proba", sep = ""))
		bad_predict = which (v_predict$class !=  Y_real)
		ACPPredict(data_apo_predict[,list_descriptor], bad_predict, Y_real, paste(path_file_descriptor, "_badPredicted"), "Bad_predicted")
	}
	else {
		print ("************************************")
		
	}
}

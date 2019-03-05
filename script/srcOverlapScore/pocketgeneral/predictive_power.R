#!/usr/bin/env Rscript

# By BORREL Alexandre
# 04-2012


# prend 2 listes et sa calcul TP, TN, FP et FN
calculTaux = function (list_predict, list_real){
	nb_value = length (list_real)
	i = 1
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	while(i <= nb_value ){
		if (as.character(list_predict[i])=="d"){
			if (list_predict[i] == list_real[i]){
				tp = tp + 1
			}else {
				fp = fp + 1
			}
		}else{
			if (list_predict[i] == list_real[i]){
				tn = tn + 1
			}else {
				fn = fn + 1
			}
		}
		i = i + 1
	}
	print (paste ("TP : ", tp, sep = ""))
	print (paste ("TN : ", tn, sep = ""))
	print (paste ("FP : ", fp, sep = ""))
	print (paste ("FN : ", fn, sep = ""))

	taux = c(tp,tn,fp,fn)
	return (taux)
}

# Calcul rate without print values
calculTaux2 = function (list_predict, list_real){
	
	nb_value = length (list_real)
	i = 1
	tp = 0
	fp = 0
	tn = 0
	fn = 0
	while(i <= nb_value ){
		if (as.character(list_predict[i])=="d"){
			if (list_predict[i] == list_real[i]){
				tp = tp + 1
			}else {
				fp = fp + 1
			}
		}else{
			if (list_predict[i] == list_real[i]){
				tn = tn + 1
			}else {
				fn = fn + 1
			}
		}
		i = i + 1
	}
	taux = c(tp,tn,fp,fn)
	return (taux)
}




######ACC PRECISION  RECALL##########

accuracy = function (tp, tn, fp, fn){
	return ((tp + tn)/(tp + fp + tn +fn))
}

precision = function (tp, fp){
	return (tp/(tp + fp))
}

recall = function (tp, fn){
	return (tp/(tp + fn))
}

specificity = function (tn, fp){
	return (tn/(tn + fp))
}

sensibility = function (tp, fn){
	return (tp/(tp + fn))
}


BCR = function (tp, tn, fp, fn){
	return (0.5*(tp/(tp+fn) + tn/(tn+fp)))

}


MCC = function (tp, tn, fp, fn){
	numerator = tp*tn-fp*fn
	denumerator = (tp+fp) * (tp+fn) * (tn+fp) * (tn+fn)
	return (numerator / sqrt(denumerator))
}


qualityPredict = function (predict, Y2){
	print (as.vector(predict)[[1]])
	print (as.vector(Y2)[[1]])
	v_predict = calculTaux (as.vector(predict)[[1]], as.vector(Y2)[[1]])
	print (paste ("accuracy : ", accuracy(v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
	print (paste ("precision : ",precision(v_predict[1], v_predict[3]), sep = ""))
	#print (paste ("recall : ", recall(v_predict[1], v_predict[4]), sep = ""))
	print (paste ("sensibility : ", sensibility(v_predict[1], v_predict[4]), sep = ""))
	print (paste ("specificity : ", sensibility(v_predict[2], v_predict[3]), sep = ""))
	print (paste ("BCR (balanced classification rate) : ", BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
	print (paste ("BER (balanced error rate) : ", 1 - BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
	print (paste ("MCC (Matthew) : ", MCC (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
	return (v_predict)
}




qualityPredictList = function (test_vector, real_vector){
	v_predict = calculTaux (test_vector, real_vector)
	print (paste ("accuracy : ", accuracy(v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
	print (paste ("precision : ",precision(v_predict[1], v_predict[3]), sep = ""))
	#print (paste ("recall : ", recall(v_predict[1], v_predict[4]), sep = ""))
	print (paste ("sensibility : ", sensibility(v_predict[1], v_predict[4]), sep = ""))
	print (paste ("specificity : ", sensibility(v_predict[2], v_predict[3]), sep = ""))
	print (paste ("BCR (balanced classification rate) : ", BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
	print (paste ("BER (balanced error rate) : ", 1 - BCR (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
	print (paste ("MCC (Matthew) : ", MCC (v_predict[1], v_predict[2], v_predict[3], v_predict[4]), sep = ""))
	return (v_predict)
}









cumulTaux = function (taux1, taux2){

	tp = taux1[1] + taux2[1]
	tn = taux1[2] + taux2[2]
	fp = taux2[3]
	fn = taux2[4]

	print (paste ("accuracy : ", accuracy(tp, tn, fp, fn), sep = ""))
	print (paste ("precision : ",precision(tp, fp), sep = ""))
	#print (paste ("recall : ", recall(tp, fn), sep = ""))
	print (paste ("sensibility : ", sensibility(tp, fn), sep = ""))
	print (paste ("specificity : ", sensibility(tn, fp), sep = ""))
	print (paste ("BCR (balanced classification rate) : ", BCR (tp, tn, fp, fn), sep = ""))	
	print (paste ("BER (balanced error rate) : ", 1 - BCR (tp, tn, fp, fn), sep = ""))
	print (paste ("MCC (Matthew) : ", MCC (tp, tn, fp, fn), sep = ""))	

}


# for ROC curve -> calcul vecteur prediction with probability (just for druggability)
generateVect = function(proba_out_predict, threshold){

	proba_class1 = proba_out_predict[,1]
	
	vect_out = NULL

	for (proba in proba_class1){
		if (proba > threshold){
			vect_out = c(vect_out, "d")
		}else{
			vect_out = c(vect_out, "nd")
		}
	}
	return (vect_out)
}


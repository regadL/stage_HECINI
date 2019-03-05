#!/usr/bin/env Rscript

# By BORREL Alexandre
# 04-2012

library (MASS)
require(plotrix)

#######################
#  GRAPHICS MANAGERS  #
#######################


factorACP = function (coor_point, vector_arrows){

	factor = 1
	orgin = vector_arrows
	while (max (vector_arrows[,1]) < max (coor_point[,1]) && max (vector_arrows[,2]) < max (coor_point[,2]) && min (vector_arrows[,1]) > min (coor_point[,1]) && min (vector_arrows[,2]) > min (coor_point[,2]) ){
		factor = factor + 1
		vector_arrows[,1] = vector_arrows[,1] + orgin[,1]
		vector_arrows[,2] = vector_arrows[,2] + orgin[,2]

	}
	return (factor-1)	

}

colorACPByTypeOfDescriptors = function (){
	d.col = read.csv ("temp_color", sep = "\t", header = TRUE)
	#print (d.col)
	return (d.col)
}



calculLimListData = function (data1, data2){

	x1 = min (data1)
	x2 = min (data2)
	X1 = max(data1)
	X2 = max (data2)
	l = c (min (x1,x2), max (X1,X2))
	return (l)
}


MDSElimcor = function (data, out_elimcor, path_file){

	groupe_elimcor = out_elimcor$groupes
	descriptor_selected = out_elimcor$possetap

	data = na.omit (data)
	MC = cor(data)
	dist1 = abs(1-MC)
	
	color_desc = rep("#000000",dim(data)[2])
	
	col_temp = sample(rainbow (dim(data)[2]))
	i_col = 1
	num_desc = rep(1,dim(data)[2])
	for (element in groupe_elimcor){
		color_desc[element] = col_temp[i_col]
		num_desc[element] = i_col
		i_col = i_col + 1
	}
	name_descriptor = colnames(data)

	
	png (paste (path_file, ".png", sep = ""), 2480, 1700)
	par( mar=c(12,12,12,12))
	fit <- cmdscale(as.dist(dist1),eig=TRUE, k=2)
	#c = as.vector (col.desc[,colnames (data)])
	plot (fit$points[,1], fit$points[,2], main="MDS Descriptor", xlab = "DIM 1", ylab = "DIM 2", cex.lab = 2.75, cex.axis = 2, cex.main = 2.75, type = "n", xlim = c(-1.3, 1.3), ylim = c(-1, 1))
	text (fit$points[,1], fit$points[,2]+0.05, labels = num_desc,  cex = 1.8, col = color_desc, font = 12)
	text (fit$points[descriptor_selected,1], fit$points[descriptor_selected,2]+0.08, labels = "*",  cex = 4, col = color_desc[descriptor_selected])
	text (fit$points[,1], fit$points[,2], labels = name_descriptor,  cex = 2.6, col = color_desc)
	
	dev.off()
}



# distribution
histData = function (data1, data2, path_pdf){

	nb_descriptor = dim (data1)[2]
	pdf (path_pdf)
	for (i_descriptor in 1:nb_descriptor){
		l = list (data1[,i_descriptor], data2[,i_descriptor])
		multhist(l, col = c(2,1), cex.names = 1.5, freq = TRUE, cex.axis = 1.5)
		title(main=colnames (data1)[i_descriptor], cex.main = 1.5, ylab = "Frequency", cex.lab = 1.5)
		legend("topright", col=c(2,1), legend = c("Class1", "Class2"), pch=c(16,16), lty=c(0,0,1,1,1), cex = 1.5)
	}
	dev.off()
}





# check color vector
CheckColorVector = function (l_descriptor, col.des){
	for (d in l_descriptor){
		if (is.integer0(which(d == names(col.des))) == TRUE){
			out = rep (1,length (l_descriptor))
			names (out) = l_descriptor	
			return (out)
		}
	}
	return (col.des)
}




####################
# COMPARISON TEST  #
####################


# choice parametric or non parametric test comparison
conditionTtest = function (data1, data2){
	if(sd(data1) == 0 || sd(data2) == 0){
		return (2)
	}
	
	if (length (data1) < 30 || length (data2) < 30 ){
		#normalisation
		pval1 = shapiro.test (data1)$p.value
		pval2 = shapiro.test (data2)$p.value
		
		if (pval1 < 0.005 || pval2 < 0.005){
			return (0)
		}	
	}
	pval_var = var.test (data1, data2)$p.value
	if (pval_var < 0.05){
		return (0)
	}
	return (1)
}

# run test comparison
comparisonTest = function (vector_value1, vector_value2, type){

	if (type == "parametric"){
		result = t.test (vector_value1, vector_value2)
		
	}else{
		result = wilcox.test (vector_value1, vector_value2)	
	
	}
	return (result$p.value)
}


# retrieve more significative descriptor

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



#################
# DATA MANAGERS #
#################

# separate data with class group
separeData = function (data, descriptor_class){

	data1 = data [which(data[,descriptor_class] == 0),]
	data2 = data [which(data[,descriptor_class] == 1),]

	return (list (data1, data2))
}


openData = function (path_file1, elimcor){
	desc = read.table (path_file1, header = TRUE, sep = "\t")
	#print (desc)

	# deleted line with NA
	#rownames (desc) = seq (1, dim(desc)[1])
	desc = na.omit(desc)
	 
	# dell when sd = 0
	
	sd_desc =apply (desc[,1:(dim(desc)[2])-1], 2, sd)
	

	#print (sd_desc)
	#print ("--------")
	sd_0 = which (sd_desc == 0)

	#print (sd_0)

	#print ("------------")
	#print (mode(sd_0))
	#print (length (sd_0))
	#print ("------------")
	if (length(sd_0) != 0){
		#print (as.factor (sd_0))
		#desc = desc[,-sd_0]
		desc=subset(desc,select=-sd_0)
		#print(dim(desc_new))	
	}
	if (elimcor != 0){
		drugg = desc$drugg
		desc = desc[,-which (colnames(desc) == "drugg")]
		out_elimcor = elimcor_avecY (desc, drugg, elimcor)
		descriptor = out_elimcor$possetap	
		
		MDSElimcor (desc, out_elimcor, paste (path_file1, "_", elimcor, sep = ""))

		descriptor = colnames (desc) [descriptor]
		desc = cbind (desc[,descriptor], drugg)
		return (list(desc,descriptor))
	}

	

	return (list((desc),colnames (desc)[-which (colnames(desc) == "drugg")]))
}


openDataDescriptor = function (path_file1, list_descriptor){

	#print (list_descriptor)
	desc = read.table (path_file1, header = TRUE, sep = "\t")
	list_descriptor = c (list_descriptor, "drugg")
	desc = desc[,list_descriptor]


	return (list((desc),colnames (desc)[-which (colnames(desc) == "drugg")]))
}



# change 0 or 1 by d and nd
changeList = function (list_element){
	Y = as.factor(list_element)
	Y2 = rep("d",length(Y))
	Y2[which(Y==0)]="nd"
	return (Y2)
}


# del colum
delCol = function (data, list_name_col){
	for (name_col in list_name_col){
		data = data[,-which(colnames(data)==name_col)]
	}
	return (data)
}


is.integer0 <- function(x)
{
  is.integer(x) && length(x) == 0L
}



#######################
# order data rownames #
#######################


orderByType = function (d.col, d){
	l_temp = NULL

	for (desc in colnames (d.col)){
		if (is.integer0 (which(rownames (d) == desc))== FALSE){
			l_temp = append (l_temp, desc)
		}
	}
	if (length (l_temp) != length (rownames (d))){
		print ("TO DO -> return list desc")
		return (l_temp)
	}
	d = d[l_temp,]
	#names (d) = l_temp
	return (d)
}

orderByTypeOneCol = function (d.col, d){
	l_temp = NULL

	for (desc in colnames (d.col)){
		if (is.integer0 (which(rownames (d) == desc))== FALSE){
			l_temp = append (l_temp, desc)
		}
	}
	if (length (l_temp) != length (rownames (d))){
		print ("TO DO -> return list desc")
		return (l_temp)
	}
	d = d[l_temp,]
	names (d) = l_temp
	return (d)
}



######################
# Normalization LDA  #
######################


normalizationCoef = function (coef, data_train){

	d_class1 = data_train[which(data_train[,"drugg"]==0),]
	d_class2 = data_train[which(data_train[,"drugg"]==1),]

	m_class1 = mean (d_class1[,1])
	m_class2 = mean (d_class2[,1])
	
	v_c1 =  sum((d_class1[,1]-m_class1) * (d_class1[,1]-m_class1))
	v_c2 =  sum((d_class2[,1]-m_class2) * (d_class2[,1]-m_class2))

	
	v_out = sqrt((v_c1 + v_c2) / (dim(data_train)[1] - 2))

	return (coef*v_out)

}


normalizationScalingLDA = function (scalingLDA, d){

	print ("****")
	print (colnames (d))
	print ("-----")

	l_out = NULL
	l_des = names (scalingLDA)
	
	for (desc in l_des){
		l_out = append (l_out, normalizationCoef (scalingLDA[desc], d[,c(desc,"drugg")]))
	}
	names (l_out) = names (scalingLDA)
	return (l_out)


}






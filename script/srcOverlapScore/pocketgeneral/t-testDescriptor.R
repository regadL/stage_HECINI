#!/usr/bin/env Rscript

# By BORREL Alexandre
# 04-2012

source ("tool.R")



########
# MAIN #
########

args <- commandArgs(TRUE)

file_drugg = args[1]
file_no_drugg = args[2]


data_drugg = read.table (file_drugg, sep = "")
data_nodrugg = read.table (file_no_drugg, sep = "")

result = as.integer (conditionTtest (as.vector (data_drugg[,1]),as.vector (data_nodrugg[,1])))

if (result == 1){
	R = t.test (data_drugg[,1], data_nodrugg[,1])
}else{
	#print ("wilcox")
	R = wilcox.test (data_drugg[,1], data_nodrugg[,1])	
}

# P value test
R$p.value

# mean
mean (data_drugg[,1])
mean (data_nodrugg[,1])

# SD
sd (data_drugg[,1])
sd (data_nodrugg[,1])


if (dim (data_drugg)[1] == dim (data_nodrugg)[1]){
	cor_linear = cor.test (data_drugg[,1], data_nodrugg[,1])
	cor = as.vector(cor_linear$estimate)
	print (cor)
	print (cor_linear$p.value)
}else{
	print(NA)
}


#!/usr/bin/env Rscript


# By BORREL Alexandre
# 11-2012

library (MASS)
source ("predictive_power.R")
source ("elimcor_avecY.R")
source ("tool.R")



# MAIN

args <- commandArgs(TRUE)
path_filin = args[1]
group_variable = args[2]

data_open = openData (path_filin, 0)
data_desc = data_open[[1]]
list_descriptor = data_open[[2]]


# remove class descriptor

for (desc in list_descriptor){
	print (paste("-------", desc, "----------", sep = " "))
	print (paste("Shapiro", shapiro.test (data_desc[,desc])$p.val, sep = " -> "))
	print (paste ("Bartlett", bartlett.test (data_desc[,desc]~data_desc[,group_variable])$p.val, sep = " -> "))
}










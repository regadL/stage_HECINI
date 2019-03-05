#!/usr/bin/env Rscript

# By BORREL Alexandre
# 11-2012

source ("tool.R")

############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_filin = args[1]
path_directory = args[2]
class = as.character(args[3])

data_global = openData (path_filin, 0)[[1]]
#print (data_global$class)
data_drugg = data_global[which(data_global[,which(colnames (data_global) == class)] == 1),]
data_nodrugg = data_global[which(data_global[,which(colnames (data_global) == class)] == 0),]

#print (data_drugg)
#print (data_nodrugg)

histData (data_drugg, data_nodrugg, paste(path_directory, "distributionDescriptor.pdf", sep = ""))




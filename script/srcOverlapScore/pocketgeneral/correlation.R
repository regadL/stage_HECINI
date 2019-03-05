#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
path_file = args[1]


data=read.table (path_file, sep = "\t", header = TRUE)

cor_linear = cor.test (data[,2], data[,3])
png (paste(path_file, ".png", sep = ""))
plot (data[,2], data[,3], main = paste ("Coef correlation\n",as.character(cor_linear$estimate), sep = ""), xlab = colnames(data)[2], ylab = colnames(data)[3], type = "n")
text (data[,2], data[,3], label = data[,1])
dev.off()



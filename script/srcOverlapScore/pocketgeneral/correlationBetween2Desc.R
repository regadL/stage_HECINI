#!/usr/bin/env Rscript
require(graphics)


args <- commandArgs(TRUE)
path_filin = args[1]

data = read.table (path_filin, sep = "\t", header = TRUE)

cor_linear = cor.test (data[,1], data[,2])


postscript (file = paste(path_filin, ".ps", sep = ""))
par(mar=c(8,8,1,8))
#plot (data[,1], data[,2], main = paste ("Coef correlation\n",as.character(cor_linear$estimate), sep = ""),xlab = colnames(data)[1], ylab = colnames(data)[2], col = data[,3] + 1, pch = 19,cex.lab = 1.8, cex.main = 2, cex = 2)

plot (data[,1], data[,2], xlab = colnames(data)[1], ylab = colnames(data)[2], col = data[,3] + 1, pch = 19,cex.lab = 1.8, cex.main = 2, cex = 2)

dev.off()


png (paste(path_filin, ".png", sep = ""),800, 800)
par(mar=c(8,8,8,8))
#plot (data[,1], data[,2], main = paste ("Coef correlation\n",as.character(cor_linear$estimate), sep = ""),xlab = colnames(data)[1], ylab = colnames(data)[2], col = data[,3] + 1, pch = 19,cex.lab = 1.8, cex.main = 2, cex = 2)

plot (data[,1], data[,2], xlab = colnames(data)[1], ylab = colnames(data)[2], col = data[,3] + 1, pch = 19,cex.lab = 1.8, cex.main = 2, cex = 2)

dev.off()


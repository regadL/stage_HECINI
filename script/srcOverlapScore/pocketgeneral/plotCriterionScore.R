#!/usr/bin/env Rscript

source ("predictive_power.R")

args <- commandArgs(TRUE)
path_file = args[1]


data = read.table (path_file, sep = "\t")


rownames (data) = data[,1]
data = data[,-1]

#print (data)

png(paste(path_file, ".png", sep = ""), 1600, 800)
par(mfrow = c(2,1))

plot (seq(1,dim(data)[1]), data[,2], type = "l", axe = FALSE, xlab = "", ylab = "", main = "MCC", ylim = c(0.5,1))
axis (2)
axis (1,seq(1,dim(data)[1]), label = rownames (data), las = 2)
points (seq(1,dim(data)[1]), data[,3], type = "l", col = "red")
segments (1:dim(data)[1],0, 1:dim(data)[1],1)
legend("topright", col=c(1,2), legend = c("LOO", "TEST"), pch=c(16,16), lty=c(0,0,1,1,1), cex = 1.5)


plot (seq(1,dim(data)[1]), data[,1], type = "l", axe = FALSE, xlab = "", ylab = "", col = "blue", main = "Criterion")
axis (2)
axis (1,seq(1,dim(data)[1]), label = rownames (data), las = 2)
segments (1:dim(data)[1],0, 1:dim(data)[1],10)

dev.off()


#!/usr/bin/env Rscript



args <- commandArgs(TRUE)
path_file = args[1]


data = read.table (path_file, sep = "\t", header = TRUE)
print (data)
rownames (data) = data[,1]
#data=data[order(data[,2],decreasing = T),]



png(paste(path_file, ".png", sep = ""), 1400, 400)


plot (seq(1,dim(data)[1]), data[,2], type = "l", axe = FALSE, xlab = "", ylab = "", col = "blue", main = "Overlap", ylim = c(0,1))
lines (seq(1,dim(data)[1]), data[,3], col = "red")
lines (seq(1,dim(data)[1]), data[,4], col = "green")
axis (2)
axis (1,seq(1,dim(data)[1]), label = rownames (data), las = 2)
segments (1:dim(data)[1],0, 1:dim(data)[1],1)

legend("topleft", col=c("blue","red","green"), legend = c("PO", "MO","RO"),pch=c(20,20,20),lty=c(1,1,1), cex = 1)
	

dev.off()


# hist commun
png(paste(path_file, "_hist.png", sep = ""), 800, 800)
par(mfrow = c(2,2))

hist(data[,5])
hist(data[,5]/data[,6])
hist(data[,5]/data[,7])

dev.off()

print ("=====")
print ("MEANS")
print (mean (data[,2]))
print (sd(data[,2]))

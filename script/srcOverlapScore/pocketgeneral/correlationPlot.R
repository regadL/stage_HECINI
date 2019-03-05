#!/usr/bin/env Rscript

# BORREL Alexandre
# 10-2013


args <- commandArgs(TRUE)

p_filin = args[1]
thresold = as.double (args[2])


d= read.table (p_filin, header = TRUE)

d_temp = d[which(d[,1] >= thresold), ]
cor.val = cor(d_temp[,2], d_temp[,3])
(p.val = cor.test (d_temp[,2], d_temp[,3])$p.value)
print (cor.test (d_temp[,2], d_temp[,3]))
print (p.val)


# png
png (paste(p_filin, ".png", sep = ""))
plot (d_temp[,2], d_temp[,3], main = paste("cor =", cor.val, "p-val =",p.val , "dim =", dim(d_temp)[1], sep = " "), xlab = colnames (d_temp)[2], ylab = colnames (d_temp)[3])
dev.off ()


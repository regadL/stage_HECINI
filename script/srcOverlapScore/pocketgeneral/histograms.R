#!/usr/bin/env Rscript

# BORREL Alexandre
# 04-2012

args <- commandArgs(TRUE)

file = args[1]
type = args[2]
brk = as.integer (args[3])
organised = as.integer (args[4])


# histograms

d = read.table (file, header = FALSE)

# cut function number col
nb_hist = dim (d)[2]


png (paste (file, ".png", sep = ""), 400, 400*nb_hist)
par (mfrow = c(nb_hist, 1))


for (i in seq (nb_hist)){
	hist (d[,i], xlim = c(min(d[,i]), max(d[,i])), breaks = brk, main = type)
}
dev.off()



# boxplot
d = d[order (d, decreasing = TRUE),]
png (paste (file, "boxplot.png", sep = ""))
if (organised == 0){
	boxplot (d, main = "boxplot global")
}else {
	par (mfrow = c(1,2))
	boxplot (d , main = "boxplot global")
	boxplot (d[1:organised,], main = paste("boxplot 1 to ",organised, sep = " "))
}

dev.off()

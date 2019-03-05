#!/usr/bin/env Rscript

args <- commandArgs(TRUE)
path_file1 = args[1]
begin = as.double(args[2])
end = as.double(args[3])
step = as.double(args[4])
type = args[5]



png(paste (path_file1, ".png", sep = ""))

for (i in seq(begin,end,step)){
	path_filin = paste (path_file1, i, sep = "")
	print (path_filin)
	d= read.table (path_filin, header = FALSE, sep = "\t")
	print (d)
	if (type == "pocket"){
		x = as.vector(d[,1])
		y = as.vector(d[,2])
		plot (x, y, col = i*100, xlim = c(5,7.0), ylim = c(10, 120), main = "Nb pocket fonction de M", xlab = "Parametre M", ylab = "Nombre de poches", type = "l")
		text (x[length(x)] + 0.1, y[length(y)] + 0.1, paste("m ",i, sep = ""))
		par(new = T)

	}
	else {
		x = as.vector(d[,1])
		y = as.vector(d[,2])
		plot (x, y, col = i*100, xlim = c(5,7.0), ylim = c(450, 1400), main = "Volume fonction de M", xlab = "Parametre M", ylab = "Volume", type = "l")
		text (x[length(x)] + 0.1, y[length(y)] + 0.1, paste("m ",i, sep = ""))
		par(new = T)
	}
}

dev.off()


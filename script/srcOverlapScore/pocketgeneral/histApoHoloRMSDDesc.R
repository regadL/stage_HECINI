#!/usr/bin/env Rscript

# BORREL Alexandre
# 01-2013


source ("tool.R")



args <- commandArgs(TRUE)

path_file_desc1 = args[1]
path_file_desc2 = args[2]
path_file_RMSD = args[3]
path_dir = args[4]
descriptor = as.character(args[5])

print (path_file_desc1)


# run
desc1 = openData (path_file_desc1, 0)[[1]]
desc2 = openData (path_file_desc2, 0)[[1]]


RMSD = read.table (path_file_RMSD, sep = "\t", header = FALSE)


print (descriptor)

print (desc1[,which(colnames(desc1)==descriptor)])

y_min = min (c (desc1[,descriptor], desc2[,descriptor],RMSD[,3]))
y_max = max (c (desc1[,descriptor], desc2[,descriptor], RMSD[,3]))


data_plot = NULL
nb_line = dim (RMSD)[1]
for (i in seq (1, nb_line)){
	new_line = c (desc1[which (rownames (desc1) == RMSD[i,1]),descriptor], desc2[which (rownames (desc2) == RMSD[i,2]),descriptor], RMSD[i,3])
	data_plot = rbind (data_plot, new_line)
}


png (paste (path_dir, descriptor, ".png", sep = ""), 1700, 1500)
par(mar=c(8,8,8,8))
plot (data_plot[,3], ylim = c(y_min,y_max), xlim = c(0,nb_line), type = "l", col ="red", xaxt = "n", main = descriptor)

points (data_plot[,1], type = "l", col = "blue")
points (data_plot[,2], type = "l", col = "orange")

grid(lwd = 1, lty = 1)
axis(1,seq(1,nb_line), RMSD[,1], cex.axis = 1.75, las = 2)
legend("topleft", col=c("red", "blue", "orange"), legend = c("RMSD", "APO","Holo") ,pch=c(20,20,20),lty=c(0,0,0,0), cex = 3.5) 


dev.off ()










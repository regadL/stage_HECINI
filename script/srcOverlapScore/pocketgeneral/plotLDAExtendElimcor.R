#!/usr/bin/env Rscript

# By BORREL Alexandre
# 04-2012

args <- commandArgs(TRUE)
path_file = args[1]
carac = as.integer(args[2])

data = read.table (path_file, sep = "\t", header = FALSE)

png (paste (path_file, ".png", sep = ""), 2480, 3508)
par( mar=c(12,12,12,12))
par(mfrow = c(2,1))

plot (data[,1],data[,2], xlab = "Correlation coefficient", ylab = "Number of Descripteurs", type = "l", cex.lab = 2.75, lwd = 3, cex.axis = 2)

if (carac == 1){
	plot (data[,1],data[,3], xlab = "", ylab = "", type = "l", xlim = c(0,1), ylim = c (0.4, 1), col = "red", lwd = 3, yaxt = "n", xaxt = "n")
	par (new = TRUE)
	plot (data[,1],data[,4], xlab = "", ylab = "", type = "l", xlim = c(0,1), ylim = c (0.4, 1), col = "chocolate1", lwd = 3, yaxt = "n", xaxt = "n")
	par (new = TRUE)
	plot (data[,1],data[,6], xlab = "", ylab = "", type = "l", xlim = c(0,1), ylim = c (0.4, 1), col = "black", lwd = 3, yaxt = "n", xaxt = "n")
	par (new = TRUE)
	plot (data[,1],data[,7], xlab = "Correlation coefficient", ylab = "Percentage", type = "l", xlim = c(0,1), ylim = c (0.4, 1), col = "blue", cex.lab = 2.75, lwd = 3, cex.axis = 2)
	legend ("topleft", legend = c("Accuracy","Precision", "Sensibility", "Specificity"),col = c("red", "chocolate1", "black", "blue" ), cex = 2.75, bty="y", lty=1)
	dev.off()
}else{

	plot (data[,1],data[,3], xlab = "", ylab = "", type = "l", xlim = c(0,1), ylim = c (0.4, 1), col = "red", lwd = 3, yaxt = "n", xaxt = "n")
	par (new = TRUE)
	plot (data[,1],data[,4], xlab = "", ylab = "", type = "l", xlim = c(0,1), ylim = c (0.4, 1), col = "chocolate1", lwd = 3, yaxt = "n", xaxt = "n")
	par (new = TRUE)
	plot (data[,1],data[,5], xlab = "Correlation coefficient", ylab = "Percentage", type = "l", xlim = c(0,1), ylim = c (0.4, 1), col = "blue", cex.lab = 2.75, lwd = 3, cex.axis = 2)
	legend ("topleft", legend = c("Accuracy LOO","Accuracy trainning", "Accuracy test"),col = c("red", "chocolate1", "blue"), cex = 2.75, bty="y", lty=1)
	dev.off()
}

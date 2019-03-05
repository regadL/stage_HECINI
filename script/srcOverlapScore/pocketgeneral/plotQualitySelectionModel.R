#!/usr/bin/env Rscript


# By BORREL Alexandre
# 10-2013



#################################
#           MAIN                #
#################################

args <- commandArgs(TRUE)
path_file = args[1]

d = read.table (path_file, sep = "\t", header = FALSE)

# histgrams
png (paste(path_file, "_hist.png", sep = ""), 800, 800)
par (mfrow = c(3,2))

# histogramm
hist (d[,1], main = "ACC LOO", col = "blue")
hist (d[,2], main = "ACC train", col = "blue")
hist (d[,3], main = "ACC test", col = "blue")
hist (d[,4], main = "MCC LOO", col = "blue")
hist (d[,5], main = "MCC train", col = "blue")
hist (d[,6], main = "MCC test", col = "blue")
dev.off ()

# hist selected
png (paste(path_file, "_hist_select.png", sep = ""), 800, 800)
par (mfrow = c(3,2))
hist (d[which(d[,1] > 0.80),1], main = "ACC LOO", breaks = 10, col = "blue")
hist (d[which(d[,2] > 0.80),2], main = "ACC train", breaks = 10, col = "blue")
hist (d[which(d[,3] > 0.80),3], main = "ACC test", breaks = 10, col = "blue")
hist (d[which(d[,4] > 0.70),4], main = "MCC LOO", breaks = 10, col = "blue")
hist (d[which(d[,5] > 0.70),5], main = "MCC train", breaks = 10, col = "blue")
hist (d[which(d[,6] > 0.70),6], main = "MCC test", breaks = 10, col = "blue")
dev.off ()


# plot
png (paste(path_file, "_plot.png", sep = ""), 800, 800)
par (mfrow = c(2,1))

# plot ACC
plot(d[order(d[,1],decreasing = F),1] ,col = "red", main = "ACC selected model", ylim = c(0.6,1), pch = 19, cex = 0.2)
par (new = TRUE)
plot (d[order(d[,1],decreasing = F),2] ,col = "blue", ylim = c(0.6,1), pch = 19, cex = 0.2)
par (new = TRUE)
plot (d[order(d[,1],decreasing = F),3] ,col = "cyan", ylim = c(0.6,1), pch = 19, cex = 0.2)


#plot MCC
plot(d[order(d[,4],decreasing = F),4] ,col = "red", main = "MCC selected model", ylim = c(0.6,1), pch = 19, cex = 0.2)
par (new = TRUE)
plot (d[order(d[,4],decreasing = F),5] ,col = "blue", ylim = c(0.6,1), pch = 19, cex = 0.2)
par (new = TRUE)
plot (d[order(d[,4],decreasing = F),6] ,col = "cyan", ylim = c(0.6,1), pch = 19, cex = 0.2)
dev.off ()


# vs
png (paste(path_file, "combi_plot.png", sep = ""), 800, 800)
par (mfrow = c(3,2))
plot(d[order(d[,4],decreasing = F),4], d[order(d[,4],decreasing = F),5] , main = "MCC LOO Vs Train", pch = 19, cex = 0.5, xlab = "MCC LOO", ylab = "MCC Train")
plot(d[order(d[,4],decreasing = F),4], d[order(d[,4],decreasing = F),6] , main = "MCC LOO Vs Test", pch = 19, cex = 0.5, xlab = "MCC LOO", ylab = "MCC Test")
plot(d[order(d[,4],decreasing = F),5], d[order(d[,4],decreasing = F),6] , main = "MCC Train Vs Test", pch = 19, cex = 0.5, xlab = "MCC Train", ylab = "MCC Test")


plot(d[order(d[,1],decreasing = F),1], d[order(d[,4],decreasing = F),2] , main = "ACC LOO Vs Train", pch = 19, cex = 0.5, xlab = "ACC LOO", ylab = "ACC Train")
plot(d[order(d[,1],decreasing = F),1], d[order(d[,4],decreasing = F),3] , main = "ACC LOO Vs Train", pch = 19, cex = 0.5, xlab = "ACC LOO", ylab = "ACC Test")
plot(d[order(d[,1],decreasing = F),2], d[order(d[,4],decreasing = F),3] , main = "ACC Train Vs Test", pch = 19, cex = 0.5, xlab = "ACC Train", ylab = "ACC Test")



dev.off ()



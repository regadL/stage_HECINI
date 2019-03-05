#!/usr/bin/env Rscript

# By BORREL Alexandre
# 05-2013

library(ggplot2)
library(plotrix)

############
#   MAIN   #
############

args <- commandArgs(TRUE)
path_file_proba_family = args[1]


d_family = read.table(path_file_proba_family, sep = "\t")
l_family = unique (d_family[,3])
l_family = l_family[order(l_family)]

col_plot = NULL
for (family in l_family ){
	i = which(d_family[,3] == family)[1]
	if(d_family[i,4] =="d"){
		col_plot = append(col_plot,2)
	}else {
		col_plot = append(col_plot,1)
	}
}

png (paste(path_file_proba_family, "_d.png"), 1200, 1200)
#par(mfrow = c(2,1))
par(mar = c(30,4,4,4))
#boxplot(d_drug[,2]~d_drug[,3], medlwd = 1, las = 2)
boxplot(d_family[,2]~d_family[,3], medlwd = 1, las = 2, col = col_plot)
abline(h=0.5)
dev.off()

png (paste(path_file_proba_family, "_nd.png"), 1200, 2400)
par(mar = c(30,4,4,4), mfrow = c(2,1))
boxplot(d_family[which(d_family[,4] == "d"),2]~d_family[which(d_family[,4] == "d"),3], medlwd = 1, las = 2, col = col_plot)
abline(h=0.5)
boxplot(d_family[which(d_family[,4] == "n"),2]~d_family[which(d_family[,4] == "n"),3], medlwd = 1, las = 2, col = col_plot)
abline(h=0.5)
dev.off()



#!/usr/bin/env Rscript

source ("predictive_power.R")

args <- commandArgs(TRUE)
path_file = args[1]


data = read.table (path_file, sep = "\t")


rownames (data) = data[,1]
data = data[,-1]
data=data[order(data[,2],decreasing = T),]

data = na.omit(data)

#print (data)
color = rep("",dim(data)[1])
#print (length(color))
#print (which (data[,2] > 0.5 & data[,4] == 1))
#print (which (data[,2] <= 0.5 & data[,4] == 2))

color [which (data[,2] > 0.5 & data[,4] == 1)] = rownames(data)[which (data[,2] > 0.5 & data[,4] == 1)]
color [which (data[,2] <= 0.5 & data[,4] == 2)] = rownames(data)[which (data[,2] <= 0.5 & data[,4] == 2)]

png(paste(path_file, ".png", sep = ""), 1600, 800)
par(mfrow = c(2,1))

plot (seq(1,dim(data)[1]), data[,2], type = "l", axe = FALSE, xlab = "", ylab = "", main = "druggability Score")
axis (2)
axis (1,seq(1,dim(data)[1]), label = rownames (data), las = 2)
axis (1,seq(1,dim(data)[1]), label = color, las = 2, col.axis = "red")
points (seq(1,dim(data)[1]), data[,3], type = "l", col = "red")
segments (1:dim(data)[1],0, 1:dim(data)[1],1)
segments (0, 0.5, dim(data)[1], 0.5)

plot (seq(1,dim(data)[1]), data[,1], type = "l", axe = FALSE, xlab = "", ylab = "", col = "blue", main = "Overlap")
axis (2)
axis (1,seq(1,dim(data)[1]), label = rownames (data), las = 2)
axis (1,seq(1,dim(data)[1]), label = color, las = 2, col.axis = "red")
segments (1:dim(data)[1],0, 1:dim(data)[1],1)

dev.off()


real_vector = data[,4]
real_vector[which (real_vector == 2)] = "d"
real_vector[which (real_vector == 1)] = "nd"

dogstie_vector = data[,3]
dogstie_vector[which (dogstie_vector >= 0.5)]  = "d"
dogstie_vector[which (dogstie_vector < 0.5)]  = "nd"

drugpred_vector = data[,2]
drugpred_vector[which (drugpred_vector >= 0.5)]  = "d"
drugpred_vector[which (drugpred_vector < 0.5)]  = "nd"

print ("##############")
print ("SCORE DRUGPRED")
print ("##############")
qualityPredictList(drugpred_vector, real_vector)
print ("##############")
print ("SCORE TESTMODEL")
print ("##############")
qualityPredictList(dogstie_vector, real_vector)




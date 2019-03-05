#!/usr/bin/env Rscript





plotScore = function (data_in, path_file_png, type){

	### workaround so that lattice does not order bank names alphabetically
	
	matrix_plot = as.matrix(rbind(data_in$ScoreWithLigand, data_in$ScoreWithoutLigand))

	color_hist = as.matrix(rbind(data_in$drugg, data_in$drugg))
	
	png(path_file_png, 2480, 3508 )
	par( mar=c(10,10,10,10))

	barplot(matrix_plot, beside = TRUE, horiz = TRUE, col = color_hist, names.arg = data_in$PDB, xlim = c(0,1), las = 2, cex.axis = 2, cex.names = 2.8)
	grid(2, NA, col = 1, lwd = 10, lty = 1) 
	legend ("topright", legend=c("With ligand", "without ligand"), cex = 2.8)

	dev.off()

}



diviseDataPredicted = function(data_in){

	data_good = NULL
	data_bad = NULL
	data_goodbad = NULL


	nb_line = dim (data_in)[1]
	for (i in seq(nb_line)){
		if (data_in[i,4] == 2){
			if (data_in[i,2] >= 0.5 && data_in[i,3] >= 0.5){
				data_good = rbind (data_good, data_in[i,])
			}else if (data_in[i,2] <= 0.5 && data_in[i,3] <= 0.5){
				data_bad = rbind (data_bad, data_in[i,])
			}else{
				data_goodbad = rbind (data_goodbad, data_in[i,])
			}
		}else{
			if (data_in[i,2] < 0.5 && data_in[i,3] < 0.5){
				data_good = rbind (data_good, data_in[i,])
			}else if (data_in[i,2] > 0.5 && data_in[i,3] > 0.5){
				data_bad = rbind (data_bad, data_in[i,])
			}else{
				data_goodbad = rbind (data_goodbad, data_in[i,])
			}

		}
	}
	return (list (data_good, data_bad, data_goodbad))
}



#############
#   MAIN    #
#############


args <- commandArgs(TRUE)
path_file = args[1]
path_dir = args[2]
type = args[3]


data_score = read.table (path_file, sep = "\t", header = FALSE)
colnames (data_score) = c("PDB", "ScoreWithLigand", "ScoreWithoutLigand", "drugg")
# 1-> non druggable
# 2-> druggable

data_divised = diviseDataPredicted (data_score )
data_good = data_divised[[1]]
data_bad = data_divised[[2]]
data_goodbad = data_divised[[3]]


plotScore (data_good, paste(path_dir, "good_predict.png", sep = ""))
plotScore (data_bad, paste(path_dir, "bad_predict.png", sep = ""))
plotScore (data_goodbad, paste(path_dir, "goddbad_predict.png", sep = ""))


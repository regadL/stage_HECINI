
library(dendextend)

dend = as.dendrogram(hc)

plot(dend)

poche_number <- rep("Other", length(rownames(Mat)))


is_x <- grepl("pocket0", rownames(Mat))
poche_number[is_x] <- "0"
is_x <- grepl("pocket1", rownames(Mat))
poche_number[is_x] <- "1"

is_x <- grepl("pocket11", rownames(Mat))
poche_number[is_x] <- "11"
is_x <- grepl("pocket10", rownames(Mat))
poche_number[is_x] <- "10"


is_x <- grepl("pocket2", rownames(Mat))
poche_number[is_x] <- "2"
is_x <- grepl("pocket3", rownames(Mat))
poche_number[is_x] <- "3"
is_x <- grepl("pocket4", rownames(Mat))
poche_number[is_x] <- "4"
is_x <- grepl("pocket5", rownames(Mat))
poche_number[is_x] <- "5"
is_x <- grepl("pocket6", rownames(Mat))
poche_number[is_x] <- "6"
is_x <- grepl("pocket7", rownames(Mat))
poche_number[is_x] <- "7"
is_x <- grepl("pocket8", rownames(Mat))
poche_number[is_x] <- "8"
is_x <- grepl("pocket9", rownames(Mat))
poche_number[is_x] <- "9"



poche_number <- factor(poche_number)

length(poche_number)

n_poche_numbers <- length(unique(poche_number))



cols_11 <- colorspace::rainbow_hcl(n_poche_numbers, c = 70, l  = 50)


col_poche_number <- cols_11[poche_number]


# extra: showing the various clusters cuts 

k234 <- cutree(dend, h = c(1.3,1.5,1.6))



#names(clusters_Desc) = gsub(".desc", replacement = "", names(clusters_Desc) )

#dt_test = merge(k234,clusters_Desc, by = 0 , all = F)


# color labels by car company:

labels_colors(dend) <- col_poche_number[order.dendrogram(dend)]

# color branches based on cutting the tree into 22 clusters:

dend <- color_branches(dend, k = 23)

### plots

par(mar = c(17,4,1,1))

plot(dend)

colored_bars(cbind(k234[,], col_poche_number[]), dend, rowLabels = c(paste0("h = ", c(1.3,1.5,1.6))))

#colored_bars( as.matrix(cbind(dt_test,col_poche_number)) ,  dend)


      
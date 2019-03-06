dt = read.table("score.csv",sep = ",", header = T)
dt = dt[,2:3]
library(stringr)
dt_new = as.data.frame(str_split_fixed(dt$poches, ";", 2))
df = cbind(dt_new,dt[,2])
colnames(df)[3] = 'scores'

################## Row and colnames #####################

vec = unique(df[,1])
vec2 = unique(df[,2])

################## CREATION DE LA MATRICE ###############

Mat = matrix(NA,length(vec),length(vec2))
rownames(Mat)= vec
colnames(Mat) = vec

############################## fucntion qui cherche le score ##################


score = function(poche1,poche2,df){
  
  if(poche1 == poche2){
    s = 1
    return(s)
  }
  
  x = df[which(df[,1]== poche1) ,]
  
  if( length(x[which(x[,2] == poche2) ,3]) == 0 ){
    x = df[which(df[,1]== poche2) ,]
    s = x[which(x[,2] == poche1) ,3]
    return(s)
  } else if ( length(x[which(x[,1] == poche2) ,3]) == 0 ){
    x = df[which(df[,1]== poche1) ,]
    s = x[which(x[,2] == poche2) ,3]
    return(s)
    
  } 
}

############################ Remplissage de la Mat ##################


for (i in 1:length(vec)) {
  for (j in 1:length(vec))
    Mat[i,j] = score(rownames(Mat)[i],colnames(Mat)[j], df)
}


######################### LET S PHEATMAP THISSSSS #######################

library(pheatmap)
library(RColorBrewer)
library(viridis)

mat_colors <- list(group = brewer.pal(3, "Set1"))

map = pheatmap(Mat,  border_color      = NA,
         #color = inferno(50),
         show_colnames     = FALSE,
         show_rownames     = FALSE,
         annotation_colors = mat_colors,
         fontsize          = 8,
         main              = "HeatMap")


#######################     clusters     ######################################

poches_clust =  cbind(Mat, cluster = cutree(map$tree_row, k = 10))

#clust = as.data.frame(poches_clust[,184])



  

library(NbClust)

as.matrix(Mat)

NbClust(as.data.frame(Mat), distance = "euclidean", method = "ward.D2")



########################### des tests ... ####################################


df["3S45_pocket0_",2]

x = df[which(df[,1]== "3S45_pocket0_") ,]

x[which(x[,2] == "3ECG_pocket5_") ,3]

y = df[which(df[,1]== "3S45_pocket0_") ,]


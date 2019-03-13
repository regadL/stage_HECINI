getwd()
dt = read.table("/home/hecini/Research/stage_HECINI/script/nscore.csv",sep = ",", header = T)
dt = dt[,2:3]
library(stringr)
dt_new = as.data.frame(str_split_fixed(dt$poches, ";", 2))
df = cbind(dt_new,dt[,2])
colnames(df)[3] = 'scores'


################## Row and colnames #####################

vec = unique(df[,1])
vec2 = unique(df[,2])
length(vec)

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
         main              = "HeatMap" )


mat.dist = as.dist(1-Mat)
map = pheatmap(Mat,
#color = inferno(50),
show_colnames     = FALSE,
show_rownames     = FALSE,
fontsize          = 8,
main              = "HeatMap",
clustering_distance_rows = mat.dist,
clustering_distance_cols = mat.dist,
clustering_method = "ward.D2")

hc = hclust(as.dist(1-Mat), method="ward.D2")

plot(hc)




### score moyen par groupe ####

#par(mfrow = c(3,2))
for (c in seq(1.1 , 3 , by = 0.1)){ #des valeurs de H differentes 
groupes=cutree(hc,h=c) # je coupe 
nbr_groupes = max(groupes[]) # la nombre de clusters générés 
vect_score = NULL
for (groupe in 1:nbr_groupes){ # pour chaque cluster 
nbr_ind = length(which(sort(groupes)[] == groupe)) # je calcule le nombre d'inv par cluster
groupe_cluster = which(sort(groupes)[] == groupe) #le nom de mes individus 
s = 0
comp = 0
moy = 0
for (ind in 1:(nbr_ind-1)){
deb = ind
for(ind2 in (deb+1):nbr_ind){
s = s + score(names(groupe_cluster[ind]),names(groupe_cluster[ind2]),df)
comp = comp+1
}
}
moy = (s/comp)
vect_score = c(vect_score,moy)
}
plot(c(1:nbr_groupes), vect_score , xlab = as.character(c))
abline(h = mean(vect_score), col = 2)
print(c)
print(sum ((vect_score - mean(vect_score))^2)/length(vect_score))
}


# couper 



taille = max(groupes=cutree(hc,h=1.2))

groupes=cutree(hc,h=1.2)

sort(groupes[])



#numbre d'individus par cluster 

for (x in 1:max(groupes[])){
print(as.character(x))
print(length(which(sort(groupes[]) == x)))

}

## récupếrer les membres des clusters

nac = names(which(sort(groupes[]) == 17))
writeLines(nac, sep = " ")





for(v in 1:taille){
  
  nac = names(which(sort(groupes[]) == v))
  print(as.character(v))
  writeLines(nac, sep = " ")
  
}


# étude de la variabilité ! 


var_data = list()

for(v in 1:taille){
  
  nac = names(which(sort(groupes[]) == v))
  print(as.character(v))
  writeLines(nac, sep = " ")
  var_data[[v]] = apply(data_new[nac,], 2, sd)
  
}



##############################
##############################

big_data = as.data.frame(do.call(rbind, var_data))

par(mfrow = c(3,2))

for(b in 1:37){
  plot(c(1:24),big_data[,b], xlab = colnames(big_data)[b])
  }


head(big_data)[,1:6]


?matrix

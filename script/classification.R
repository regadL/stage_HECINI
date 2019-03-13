data =  read.table("/home/hecini/Research/stage_HECINI/data/PR2/matrice_descripteurs.mat")

which(is.na(data) == T) # pas de valeurs manquantes 

#visualisation des données  : 

boxplot(scale(data))

# nettoyage de données 

# 1/ variables avec une variance proche de zero : 

# 1ere méthode : j'ai pas toute les variables .. 

var = apply(data, 2, sd)
which(var == 0)

# 2eme méthode : PARFAIT 

library(caret)

vzero = nearZeroVar(data) # vecteur de variables avec une variance proche de zero 
data_new = data[,-vzero] # un nouveau data frame sans ces variables 
dim(data_new)

#matrice de corrélation

# créer la matrice de corrélation : 

matcor = cor(data_new)
head(matcor)

# eliminer les variables avec une correlation importante : 

vcorr = findCorrelation(matcor, cutoff = 0.75)
vcorr
data_new = data_new[,-vcorr]
dim(data_new)


# visualisation : 

library(corrplot)
dt_cor = corrplot(matcor, method = "circle", type = 'upper')


### certaines variables n ont pas été éliminé du coup je repete le nettoyage 

# eliminer les variables avec une correlation importante : 

vcorr = findCorrelation(matcor, cutoff = 0.85)
vcorr
data_new = data_new[,-vcorr]
dim(data_new)

# créer la matrice de corrélation : 

matcor = cor(data_new)

# visualisation : 

dt_cor = corrplot(matcor, method = "circle", type = 'upper')

### classification : 

## le nombre optimale de cluster?



library(factoextra)

# Elbow method : propose 3 ou 4 
fviz_nbclust(data_new, kmeans, method = "wss")

# Silhouette method : propose 4 

fviz_nbclust(data_new, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")

library(cluster)

## la méthode de GAP :  1 cluster?? 

gap_stat= clusGap(data_new, FUN = kmeans, nstart = 50,
                    K.max = 10, B = 200)
plot(gap_stat)

#####  le nombre optimal : manuellement ###########

#Elbow Method 

set.seed(123)

k.max <- 15

scaled_data = scale(data_new)

wss <- sapply(1:k.max, 
              function(k){kmeans(scaled_data, k, nstart=50,iter.max = 100 )$tot.withinss})
wss

plot(1:k.max, wss,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")


for (i in 1: length(wss)-1){
  print(wss[i]-wss[i+1])

}



###########################################



## lancer la kmeans avec  k = 3
"-------------------------------"

clusters = kmeans(scaled_data, centers = 3, nstart=100,iter.max = 500)

#HCLUST
"-----"

dt_dist = dist(scaled_data)

clah = hclust(dt_dist, method = "ward.D2")

plot(clah)

clusters_tree = cutree(clah, k = 3)


############################################

#PCA

library(factoextra)
library(FactoMineR)
acp = PCA(data_new, graph = F)

fviz_cluster(clusters, scaled_data, geom = "point")

fviz_pca_var(acp)

fviz_ca_biplot(acp)

fviz_pca_biplot(acp, 
                # Individus
                geom.ind = "point",
                fill.ind = clusters$cluster, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",
                # Variables
                select.var = list(contri=15)

)




# Compute hierarchical clustering and cut into 4 clusters
res <- hcut(scaled_data, k = 3)
# Visualize
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))

?hcut

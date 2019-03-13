data =  read.table("/home/hecini/Research/stage_HECINI/data/PR2/mat.mat")

rownames(data) = gsub(".desc", replacement = "", rownames(data) )

which(is.na(data) == T) # pas de valeurs manquantes 

dim(data)

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


#####  le nombre optimal : manuellement ###########

#Elbow Method 

set.seed(123)

k.max <- 15

sc_data = scale(data_new)

wss <- sapply(1:k.max, 
              function(k){kmeans(sc_data, k, nstart=50,iter.max = 100 )$tot.withinss})
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

clusters_k = kmeans(sc_data, centers = 3, nstart=100,iter.max = 500)

#HCLUST
"-----"

dt_distDesc = dist(sc_data)

clah_Des = hclust(dt_distDesc, method = "ward.D2")

plot(clah_Des)
rect.hclust(clah_Des, 22)
clusters_Desc = cutree(clah_Des, k = 22)


############################################

#PCA

acp_desc = PCA(data_new, graph = F)

fviz_cluster(clusters_k, sc_data, geom = "point", select.var = list(contrib = 15) )

fviz_pca_var(acp_desc)

fviz_pca_biplot(acp_desc, 
                # Individus
                geom.ind = "point",
                fill.ind = clusters_k$cluster, col.ind = "black",
                pointshape = 21, pointsize = 2,
                palette = "jco",

                # Variables
            gradient.cols = "RdYlBu",
            select.var = list(contrib = 15)
                
)


# Contribution totale sur PC1 et PC2

fviz_contrib(acp_desc, choice = "var", axes = 1:2)


#préparation des données : 

dt = read.table("barycentre_new",sep = ",", header = T)

dt = dt[,c(2,4:6)] 

dim(dt)

vec = dt[dt[,1] == "1HSI_superpos.pdb",2:4]  # x y et z de 1HSI_pocket
rownames(vec)

#calcul des distannce 

val = NULL
racine =  NULL
for (i in 1:nrow(dt)) {
  
  racine = sqrt( (dt[i,2]-vec[1])^2 + (dt[i,3]-vec[2])^2 + (dt[i,4]-vec[3])^2  )
  val = c(val,racine)
  
  
}


dt$distance = val # ajout de la nouvelle colonne 

ligne  = as.numeric(rownames(vec))
dt = dt[-ligne,] # supprimer 1HSI


# classification : 

library(factoextra)

# trouver le nombre optimal de cluster 

# Elbow method :

fviz_nbclust(dt[,-1], kmeans, method = "wss", k.max = 20) # pour moi k = 3 // 10 est interessant aussi 

# Méthode Nbclust : elle propose 3 

library(NbClust)
NbClust(data =as.data.frame(lapply(dt[,-1], unlist)) , distance = "euclidean",
        min.nc = 2, max.nc = 15, method = "kmeans")

## lancer la kmeans avec  k = 3

"-------------------------------"

dtt = as.data.frame(lapply(dt[,-1], unlist))
rownames(dtt) <- dt[,1]

clusters_pos = kmeans(dtt[,-4], centers = 3, nstart=100,iter.max = 500)

clusters_pos$cluster

#HCLUST
"-----"

dt_dist = dist(dtt[,-4])

clah = hclust(dt_dist, method = "ward.D2")

plot(clah)

clusters_tree = cutree(clah, k = 3)


############################################

#PCA
library(factoextra)
library(FactoMineR)

acp_pos = PCA(dtt[,-4], graph = F)

#visualiser les clusters :
fviz_cluster(clusters_pos, dtt, geom = "point")
fviz_pca_ind(acp_pos, col.ind = clusters_pos$cluster)

# les variables et les individus : 
fviz_pca_var(acp_pos)
fviz_pca(acp_pos , col.ind = clusters_pos$cluster )

c1 = which(clusters_pos$cluster == 1)

write.table(c1, file = "lala", append = FALSE, sep = " ", dec = ".",
            row.names = T, col.names = F)

fviz_pca_biplot(acp_pos, 
                # Individus
                geom.ind = "point",
                col.ind = clusters_pos$cluster,
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                pointsize = 2
)


cmd.res <- cmdscale(dist(dtt[,-4]))
plot(cmd.res[,1:2], col = clusters_pos$cluster, pch=19)

# Compute hierarchical clustering and cut into 4 clusters
res <- hcut(dtt, k = 3)
# Visualize
fviz_dend(res, rect = TRUE, cex = 0.5,
          k_colors = c("#00AFBB","#2E9FDF", "#E7B800"))



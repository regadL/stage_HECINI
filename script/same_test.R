clust_type1 = clusters_k$cluster
clust_type2 = clusters$cluster

t1 = as.data.frame(clust_type1)
t2 = as.data.frame(clust_type2)

rownames(t2) = paste(rownames(t2),'.desc',sep = "")
de <- merge(t1, t2, by=0, all=TRUE)


library('clusteval')
res = cluster_similarity(de$clust_type1, de$clust_type2, 
                         similarity="jaccard", method="independence")
res


#####################################################

nom = NULL
num_cluster = NULL
for (i in 1:nrow(data_comparaison)){
  if (data_comparaison[i,1] == data_comparaison[i,2]){
    nom = c(nom,rownames(data_comparaison)[i])
    num_cluster = c(num_cluster,data_comparaison[i,1])
  }
}
nom
num_cluster

same_cluster = data.frame(nom,num_cluster)

#####################################################



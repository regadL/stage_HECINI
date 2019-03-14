#calcul des scores entre les groupes : 


score_all = NULL


for (c in seq(1.1 , 3 , by = 0.1)) { #des valeurs de H differentes 
  
  groupes=cutree(hc,h=c) # je coupe 
  
  nbr_groupes = max(groupes[]) # le nombre de clusters générés 
  
  vect_score = NULL
  
  #créer la matrice 
  matrice = matrix(NA,nbr_groupes,nbr_groupes)
  rownames(matrice) = c(1:nbr_groupes)
  colnames(matrice) = c(1:nbr_groupes)

  for (groupe in 1:nbr_groupes){ # pour chaque cluster 
    
    nbr_ind = length(which(sort(groupes)[] == groupe)) #nbr de poches dans le clus 1 
    
    nom_cluster1 = which(sort(groupes)[] == groupe) #le nom des poches1

    for(groupe2 in 1:nbr_groupes){

      s = 0 #la somme
      
      comp = 0 #compteur 
      
      moy = 0 #moyenne 
      
      nbr_ind2 = length(which(sort(groupes)[] == groupe2)) #nombre de poches dans clus 2 
      
      nom_cluster2 = which(sort(groupes)[] == groupe2) #le nom des poches2
      
      for (ind in 1:nbr_ind){ # pour chaque poche1
        
        for(ind2 in 1:nbr_ind2){ # pour chaque poche2
          
          s = s + score(names(nom_cluster1[ind]),names(nom_cluster2[ind2]),df)
          
          comp = comp+1
          
        }
      }
      
      #remplir la matrice 
      
      moy = (s/comp)
      matrice[groupe,groupe2] = round(moy*100,2)
    }

  }
  
  form = sprintf('%s_clusters.csv', c)
  write.csv(matrice, file = form, row.names = F , col.names = F)
  

  }
    


#calcul de la somme d'une demi matrice : 
vec_somme = NULL

H_couper = NULL

vec_somme2 = NULL

diag_somme = NULL

diag_mean = NULL

for (z in seq(1.1 , 3 , by = 0.1)){
  
  doc = sprintf("%s_clusters.csv",z)
  
  m1 = read.csv(doc)
  
  H_couper = c(H_couper,z)
  
  somme = 0
  
  somme2 = 0
  
  diag_somme = c(diag_somme,sum(diag(as.matrix(m1))))
  diag_mean = c(diag_mean, mean(diag(as.matrix(m1))))
  
  for(q in 2:nrow(m1)){

    for (d in 1:(q-1)){
      
      somme = somme + m1[q,d]
      if(m1[q,d]>10)
        somme2=somme2+m1[q,d]
      
      }
}
    
  vec_somme2 = c(vec_somme2, somme2) 
  vec_somme = c(vec_somme, somme)

  
  }


plot(H_couper, diag_somme, pch = 18, col = 5, type = "b", xlab = "valeurs de H", ylab = "la somme")
points(H_couper, vec_somme, pch = 19, col = 6, type = "b", xlab = "valeurs de H", ylab = "la somme")
points(H_couper, vec_somme2, pch = 18, col = 3, type = "b")
points(H_couper, diag_mean, pch = 18, col = 4, type = "b")
legend("topright", 95, legend=c("Somme de la diag", "somme de la semi matrice", 
                       "somme des valeurs les plus importante de la semi matrice ", "moyenne de la diag "),
       col=c(5,6,3,4), lty=1:2, cex=0.8)

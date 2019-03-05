# This is my README
8-10-12
- Ajout d'un main pour lancer les descripteurs à partir d'un fichier poste seulement
- modification dans tool, fonctions qui ajoute des HEADER et END
15-10-12
- ajouter les fonctions qui permettent de faire une sélection de variable en générant toute
les combinaisons de 4 descripteurs
- modification des fonctions d'écrire, pb relatif à la clef data qui provient
d'une erreur à la version de modification des écriture systeme.
- modification de la fonction R type qui permet de faire plus globalement les 
barplot, gestion du cas ou il y a seulement 2 colonnes ainsi que du cas avec les 
couleur précisé en colonne 3 (reliquat du projet salt bridges)
16-10-12
- modification mineur de lancement push pour rendre disponible la dernière version du code
22-10-12
- modification des fonction qui permettent de calculer la corrélation entre les poches estimées
avec ou sans ligand, bug non géré quand ajout de l'entrée data dans les types de poches pour 
faire une étude séparer entre les poches prisent dans les deux jeux de données (Test et Train)
25-10-12
- modification du lancement du CART pour avoir le même fichier de résultats que
les autres méthodes
31-10-12
-Changement de l'ouverture des fichiers sous R -> pb relatif à la suppresion des NA 
et des variances nulles
-Barplot des variable du LDA centrée et réduite
- Ajout de la representation du plot des variable avec la fonction de separation
entre duggable et non druggable
8-11-12
- Recode SVM, modification des outils 
14-11-12
- ajout d'une fonction de calcul du ligand pour voir si violation de Lipinski
- modification sur les SVM -> ok vis a vis de la selection de variables
15-11-12
- ajout des script R pour vérifier les conditions de validité du LDA
- ajout scripts pour faire la distribution entre bien et mal prédits
- modification majeur sur le SVM
10-12-13
- modification du calcul entre les formes apo et holo
- modification du code sur la superposition entre apo et holo
- changement des ACP pour avoir les ps avec les png pour les articles et le poster
16-01-13
- ajout d'un calcul de score de recouvrement entre les poches prédites par deux approches 
15-03-13
-finalisation de l'article, ajout de parametre dans les ACP comme les males prédites
- Chagement dans la prise en charge des ouverture de fichier, bcp moins contraignant
- Calcul de l'exposition au solvant de la poche
- Calcul des rêgles de Lipinski
- Mis à jour du nom des descripteurs (proportion) et passage du comptage dans les géométrique
03-05-13
- update finlande
- update pour le jeu DD pour faire une comparaison avec Volkamer
06-05-13
- comparaison avec volkamer ok pour le parsing des fichiers serveur, html + fichiers poches
10-06-13
- reprise la fonction d'overlap des poches pour prendre les poches et comparer les estimateurs
- finition de la fonction de comparaison
18-09-13
- ajout de nouvel estimateur cavitator mais estimateur ne fonctionne pas corectement (pb sur fortran)
- ajout de dogsite en estimateur
25-10-13
- selection de varaible basee sur le MCC
- ajout mainGeo pour faire une analyse des descripteurs associé a l'enveloppe convexe
- modification des fonctions d'ouverture de R et couleur
- integration de dogsite dans les scripts
- nouveau graphe pour la selections des descripteurs

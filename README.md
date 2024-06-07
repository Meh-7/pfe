to get folder structure:
git ls-tree -r HEAD > folder_structure.txt
02/06/2024
- on doit faire varier les n (longeurs des series)
- utiliser des modeles avec matrices non diagonales
- nb d'ob est proportionnel au nombre de parametres (10 ob par param)
- 110 ob pour nous est une petite taille si on prend tous nos 11 parametres
- on fair varier les observations pour 200, 500 , 1000 par exemple
	- faire des replications pour chaque taille
	- moyenne empirique qui nous donne une idee sur le biais
	- faire des variances, nous parle d'efficacité
	- variabilité des estimateurs: erreures MADE et MSE(vitesse de convergence et consistance), MSE nous permet de comparer aussi entre les estimateurs (OLS vs MLE par ex)
- proposition (ajouter entre parentheses "estimateur ...")
	- soit psi "Le parametre du modele donné par (2.8)"... (arranger les transposées de psi dans la proposition)
	- "l'estimateur de psi par la methodes de Moindres Carrés Ordinaires  $\hat{\Psi}_{MCO}$  est donnée par les deux etapes suivantes:"
		- ecrire "etape 1"
		- dans l'etape 2 " $\hat{\Psi}_{MCO}$ qui est un MCG, est obtenu bloc par bloc grace au theoreme de Zellner (année) "
- H01-H03 $\rightarrow$ H1 et H3 (mm chose pour les lemmes)    
-  bien introduire la section d'estimation, qu'est-ce qu'on va faire, pourqoi etc
	- ou est-ce qu'on a codé etc, les logiciels utilisés etc
- utiliser les box plot avec les tableaux pour visualiser l'information

06/06/2024
- garder a la fin du memoir dans les perspective le cas multivarie et trivarie plus general, les difficulté et le obstacles qui seront rencontrés (peut etre dans les appendices ajouter des simulations et des modeles qui existent?)
- conclusion dans rapport: 
	- enlever le "predefini" dans le threshold
	- ajouter le point de rupture (deux extensions total)
	- ajouter le cas multivarié
	- decomposer conclusion(chapitre 3) sous sections bien distingués
		- etat d'avancement
		- perspectives
- notations:
	- stat: n = taille echantillon, N = replication <- utiliser ca
	- econometrie: T = taille echantillon
- il faut savoir concentrer bien sur le pourquoi: pourquoi un certain modele, pourquoi le modele de cui et zhu au lieu de liu, le delta et son importance dans la correlation negative
- pour commission:
	- rappel rapide du provleme
	- etale par chapitre 
		- chapitre 1 : voila ce qu'on a fait ((2, 3 diapos))
			- peut etre ajouter un diagramme pour la definition des lois de poisson multivaries (il faut l'ajouter dans le memoir)
		- chapitre 2: 
			- voila notre modele, ce qui a ete fait, (cui, liu)
		- chapitre 3:
			- simulation
- figures dans la simulation: caption Y_1(t) , c'est des series chrono, no pas oublier le temps

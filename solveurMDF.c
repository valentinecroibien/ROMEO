#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Définition du type Vecteur */
typedef struct 
{ int dim; /* Nombre de coefficients */
  double* tab; /* Coefficients du vecteur */
} Vecteur;

/* Fonction qui retourne un vecteur alloué de sa dimension n */
Vecteur *alloc_Vecteur(int n){
	Vecteur * V = NULL;
	assert (n>0);

	V=(Vecteur *) malloc(sizeof(Vecteur));

	V -> dim = n;
	V -> tab = (double * ) calloc(n,sizeof(double));

	if(V -> tab==NULL){printf("Plus de place \n");exit(1);}

	return V;
}

/* Procédure qui libère la mémoire allouée à un objet de type Vecteur */
void free_Vecteur(Vecteur * V){
	free(V -> tab);
	free(V);
}

/* Procédure qui affiche un objet de type Vecteur précédé d'un message */
void affiche_Vecteur(char message[256],Vecteur * V){
	int i;

	printf("%s", message);

	for(i = 0; i < V -> dim; i++){
		printf("%f \n",V -> tab[i]);
	}
}

/* Fonction qui récupère le coefficient i d'un Vecteur - getter */
double coeff_Vecteur(Vecteur *V, int i){
  return V -> tab[i];
}

/* Procédure écrivant dans un fichier un Vecteur */
void ecrit_Vecteur(Vecteur * V, char fichier[256]){
	int i;
	FILE *f;
	
	f = fopen(fichier,"a+");
	for(i = 0; i < V -> dim; i++){
		fprintf(f,"%.16f \n",V -> tab[i]);
	}
	fclose(f);
}

/* Procédure qui affecte un coefficient à un objet de type Vecteur - SETTER */
void affecte_Vecteur(double valeur, Vecteur *V, int i){
	V -> tab[i] = valeur;
}

/* Procédure réalisant la factorisation LU d'une matrice A adaptée à un système tridiagonal */
void facto_LU(Vecteur *a, Vecteur *b, Vecteur *c, Vecteur *d, Vecteur *e, Vecteur *f, int *comp){
	int i;
	int flop=0; /*Complexité de la factorisation LU */
	int n = a -> dim;

	for(i=0; i < n-1; i++ ){
		f -> tab[i] = (c -> tab[i])/(a-> tab[i]);
		flop += 1;
		d -> tab[i] = a -> tab[i];
		e -> tab[i] = b -> tab[i];
		a -> tab[i] -= d -> tab[i];
		flop += 1;
		a -> tab[i+1] -= (f -> tab[i])*(e -> tab[i]);
		flop += 2;
		b -> tab[i] -= (f -> tab[i])*(d -> tab[i]);
		flop += 2;
		c -> tab[i] -= (e -> tab[i]);
		flop += 1;
	}
	d -> tab[n-1] = a -> tab[n-1];
	
	*comp +=flop;

}

/* Fonction descente adaptée à la résolution d'un système tridiagonal */
void descente(Vecteur * f, Vecteur* B, Vecteur *x, int *comp){
	int n = B -> dim;
	int i;
	int flop = 0; /* Complexité de la descente */

	(x -> tab[0]) = (B -> tab[0]);
	for(i=1; i<n; i++){
		x -> tab[i]= (B -> tab[i]) - (f -> tab[i-1])*(x -> tab[i-1]);
		flop += 2;
	}
	*comp += flop;
}

/* Fonction remontee adaptée à la résolution d'un système tridiagonal */
void remontee(Vecteur * d, Vecteur* e, Vecteur* b, Vecteur* x, int *comp){
	int n = b -> dim;
	int i;
	int flop=0; /* Complexité de la remontée */

	x -> tab[n-1] = (b -> tab[n-1])/(d -> tab[n-1]);
	flop += 1;
	for (i=n-2; i >= 0; i--){
		x -> tab[i] = ((b -> tab[i]) - (e -> tab[i])*(x -> tab[i+1]))/(d -> tab[i]);
		flop += 3;
	}
	*comp += flop;
}

/* Fonction solveur_LU adpatée à la résolution d'un système tridiagonal */
void solveur_LU_TriDiagonal(Vecteur *a, Vecteur *b, Vecteur *c, Vecteur *B, Vecteur *x, int *comp){
	Vecteur *d = NULL;
	Vecteur *e = NULL;
	Vecteur *f = NULL;
	Vecteur *y = NULL;
	int n = a -> dim;

	d = alloc_Vecteur(n);
	e = alloc_Vecteur(n-1);
	f = alloc_Vecteur(n-1);
	y = alloc_Vecteur(n);

	facto_LU(a,b,c,d,e,f,comp);
	descente(f,B,y,comp);
	remontee(d,e,y,x,comp);
	
	free_Vecteur(d);
  	d=NULL;
  	free_Vecteur(e);
  	e=NULL;
  	free_Vecteur(f);
  	f=NULL;
}

double f(double x){
	return 0;
}

double u(double x, double epsilon){
	return (1/(1-exp(-2/sqrt(epsilon)))) * exp(-x/sqrt(epsilon)) + (1/(1-exp(2/sqrt(epsilon)))) * exp(x/sqrt(epsilon));
}

/* Fonction qui calcule la norme infinie d'un vecteur */
double norm_inf(Vecteur *x){
  double norm=0.0;
  int i;
	
  for (i=0;i<x->dim;i++){
	if (fabs(x->tab[i]) > norm) {
		norm = fabs(x->tab[i]);
	}
  }
  
  return norm;
}

/* Fonction qui calcule la norme 2 d'un vecteur */
double norm_2(Vecteur *x){
  double norm=0.0;
  int i;
    
  for (i=0;i<x->dim;i++){
    norm += pow(fabs(x->tab[i]), 2);
  }
  
  return sqrt(norm);
}

/* Fonction soustraction C = A - B */
void soustraction(Vecteur *x, Vecteur *y, Vecteur *z){
  int i;
  assert(x->dim == y->dim);
  assert(y->dim == z->dim);
  
  for(i=0; i<=x->dim; i++){
  	z->tab[i] = x->tab[i] - y->tab[i];
  }
}

/* Fonction solveur MDF*/
double Solveur_MDF_1D(double L, int N, double a, double b, double epsilon){
	int i;
	double alpha;
	double beta;
	int comp=0;
	double einf=0.0, h;
	FILE *fichier;

	/* Vecteurs définissant la matrice A */
	Vecteur *c = NULL;
	Vecteur *d = NULL;
	Vecteur *e = NULL;
	/* Vecteur Uh */
	Vecteur *Uh = NULL;
	/* Vecteur B */
	Vecteur *B = NULL;
	/* Vecteur solution exacte */
	Vecteur *U = NULL;
	/* Vecteur erreur */
	Vecteur *EINF = NULL;
	
	/* Vecteur x */
	Vecteur *X = NULL;

	c = alloc_Vecteur(N+1);
	d = alloc_Vecteur(N);
	e = alloc_Vecteur(N);
	Uh = alloc_Vecteur(N+1);
	B = alloc_Vecteur(N+1);
	U = alloc_Vecteur(N+1);
	EINF = alloc_Vecteur(N+1);
	X = alloc_Vecteur(N+1);

	/* Valeurs définissant les matrices */
	h = (L/N);
	alpha = epsilon/(h*h);
	beta = 2*alpha + 1;

	/* Assemblage des vecteurs */
	affecte_Vecteur(1,c,0);
	affecte_Vecteur(1,c,N);
	for(i=1;i<N;i++){affecte_Vecteur(beta,c,i);}
	
	affecte_Vecteur(0,d,0);
	for(i=1;i<N-1;i++){affecte_Vecteur(-alpha,d,i);}
	
	for(i=0;i<N-1;i++){affecte_Vecteur(-alpha,e,i);}
	affecte_Vecteur(0,e,N-1);
	
	for(i=0;i<N+1;i++){
		affecte_Vecteur(i*h,X,i);
	}
	ecrit_Vecteur(X,"./X.txt");
	
		
	affecte_Vecteur(a,B,0);
	affecte_Vecteur(b,B,N);
	for(i=1;i<N;i++){affecte_Vecteur(f(X -> tab [i]), B, i);}
	
	
	/* Solution exacte */
	for(i=0; i<N+1; i++){
		affecte_Vecteur(u(X -> tab [i],epsilon), U,i);
	}
	ecrit_Vecteur(U,"./U.txt");
	
	solveur_LU_TriDiagonal(c, d, e, B, Uh, &comp);
	ecrit_Vecteur(Uh,"./Uh.txt");
	
	
	soustraction(U, Uh, EINF);
	einf = norm_inf(EINF);
	
    printf("Nombre N d'intervalles :   %8.d, h = %.6f, ||Uh||_2 = %.16e, ||U-Uh||_inf = %.16e\n",N,h,norm_2(Uh),einf);
    printf("\n");
	fichier=fopen("./Convergence.txt", "a+");
 	fprintf(fichier,"%4.0d,%.8f,%e\n", N, h, einf);
 	fclose(fichier);
 	
	return einf;
}

int main(){
    int L = 1;
    double a = 1.0, b = 0.0, epsilon = 0.00001;
    Vecteur *einf=NULL;
    Vecteur *p=NULL;
 
    FILE *f;
    int m; /* Nombre d'intervalles lu */
    int i=0; /* Nombre de Nombre d'intervalles */
    int nbl=0;
    Vecteur *NI=NULL;
    
     
    f=fopen("./NI.txt", "rt"); /* Ouverture */
    if(f==NULL) { printf("Fichier inexistant ! \n"); exit(1);}
    while(1){
    nbl=fscanf(f,"%d", &m);	/* :lecture(s) */
    if (nbl==EOF){break;}
        i += 1;
    }
    fclose(f);
    NI = alloc_Vecteur(i);
    einf = alloc_Vecteur(i);
    p = alloc_Vecteur(i-1);
    i=0;
    f=fopen("./NI.txt", "rt"); /* Ouverture */
    if(f==NULL) { printf("Fichier inexistant ! \n"); exit(1);}
    while(1){
        nbl=fscanf(f,"%d", &m);	/* :lecture(s) */
        if (nbl==EOF){break;}
        affecte_Vecteur(m,NI,i);
        i += 1;
    }
    fclose(f);
        
    printf("RÉSOLUTION PAR DIFFÉRENCES FINIES DU PROBLÈME ELLIPTIQUE 1D : \n");
    printf("{-eps u''(x) + u(x) = 0 sur [0.0; 1.0], avec eps = 1.0e-05 \n");
    printf("{            u(0.0) = 1.0\n");
    printf("{            u(1.0) = 0.0\n\n");
    printf("DE SOLUTION ANALYTIQUE : \n");
    printf(" u(x) = (1/(1-exp(-2/sqrt(epsilon)))) * exp(-x/sqrt(epsilon))\n");
    printf("      + (1/(1-exp(+2/sqrt(epsilon)))) * exp(+x/sqrt(epsilon))\n\n");
    
    /* On efface le contenu des fichiers sur lesquels nous allons réécrire les nouvelles valeurs*/
    f=fopen("./Convergence.txt", "w");
    fclose(f);
    f=fopen("./U.txt", "w");
    fclose(f);
    f=fopen("./Uh.txt", "w");
    fclose(f);
    f=fopen("./X.txt", "w");
    fclose(f);
    for(i=0; i<NI->dim; i++){
        affecte_Vecteur(Solveur_MDF_1D(L, NI -> tab[i], a, b, epsilon),einf,i);
    }

    /* Calcul et affichage de l'ordre de convergence p */
    for(i=0; i<einf->dim-1; i++){
        affecte_Vecteur((1/log(2))* log(coeff_Vecteur(einf,i)/coeff_Vecteur(einf,i+1)),p,i);
        printf("%6.d vs %6.d  ---> p = %+.16e\n", (int)(coeff_Vecteur(NI,i+1)),(int)(coeff_Vecteur(NI,i)), coeff_Vecteur(p,i));
    }

    return 0;
}

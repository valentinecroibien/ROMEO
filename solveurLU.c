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
	
	affiche_Vecteur("\nd = \n", d);
	affiche_Vecteur("\ne = \n", e);
	affiche_Vecteur("\nf = \n", f);
	printf("\nComplexité algorithmique de la factorisation LU : %d \n",flop);
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
	printf("Complexité algorithmique de la descente : %d \n",flop);
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
	printf("Complexité algorithmique de la remontée : %d \n",flop);
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

int main()
{	Vecteur *a = NULL;
	Vecteur *b = NULL;
	Vecteur *c = NULL;
	Vecteur *B = NULL;
	Vecteur *x = NULL;
	int n;
	int i;
	int complexite=0;
	double coeff;

	/* Allocation des vecteurs */
	printf("\nDonnez la taille de la matrice A : n = ? \n");
	scanf("%d", &n);
	a = alloc_Vecteur(n);
	printf("\nDonnez les coefficients de a : \n");
	for(i=0; i<n; i++){
		scanf("%lf", &coeff);
		affecte_Vecteur(coeff, a,i);
	}
	b = alloc_Vecteur( n-1);
	printf("\nDonnez les coefficients de b : \n");
	for(i=0; i<n-1; i++){
		scanf("%lf", &coeff);
		affecte_Vecteur(coeff, b,i);
	}
	c = alloc_Vecteur(n-1);
	printf("\nDonnez les coefficients de c : \n");
	for(i=0; i<n-1; i++){
		scanf("%lf", &coeff);
		affecte_Vecteur(coeff, c,i);
	}
	B = alloc_Vecteur(n);
	printf("\nDonnez les coefficients de B : \n");
	for(i=0; i<n; i++){
		scanf("%lf", &coeff);
		affecte_Vecteur(coeff, B,i);
	}
	x = alloc_Vecteur(n);


	solveur_LU_TriDiagonal(a, b, c, B, x, &complexite);
	affiche_Vecteur("\nSolution du problème tridiagonal : \n", x);
	printf("\nComplexité algorithmique totale : %d\n", complexite);

	free_Vecteur(a);
  	a=NULL;
  	free_Vecteur(b);
  	b=NULL;
  	free_Vecteur(c);
  	c=NULL;
  	free_Vecteur(B);
  	B=NULL;
  	free_Vecteur(x);
  	x=NULL;
  	
	return 0;
}

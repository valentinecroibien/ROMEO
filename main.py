import os
import sys

from datetime import datetime, time
import time as ti

from numpy import log
from pylab import *

# PRÉ TRAITEMENT
print("gcc -ansi -pedantic solveurMDF.c -lm \n")

os.system("gcc -ansi -pedantic solveurMDF.c -lm")

NI = [200, 400, 800, 1600, 3200, 6400, 128000 ]
print("Nombres d'intervalles par défaut : NI = ", str(NI),"\n")

rep = input("Est-ce que la liste vous convient ? (O/N) ")
if (rep != "O" and rep != "N") :
    print("ERREUR !")
    sys.exit(0)
    
if (rep == "N") :
	rep = input("\nEntrer un entier pour la nouvelle liste : ")
	NI = []
	while (rep != "-1") :
		NI.append(int(rep))
		rep = input("Entrer un entier pour la nouvelle liste (-1 pour quitter) : ")

print("\nNouvelle liste : NI = " + str(NI) +"\n")

ptrf=open("NI.txt","w")
for i in range(0,len(NI)) :
	ptrf.write(str(NI[i]) + "\n")
ptrf.close()


# EXÉCUTION CODE C
print("Exécution du Code C ... Running \n")
os.system("./a.out")
ti.sleep(2)
print("\nExécution du Code C ... Terminée \n")




# POST TRAITEMENT
print("Lecture des fichiers provenant du code C\n")





# Lecture du fichier X.txt
fichier = open("X.txt","r")
X = []
lignes = fichier.readlines()
for ligne in lignes :
	ligne = ligne.strip()
	xi = [float(elt) for elt in ligne.split("\n")]
	X.append(xi[0])
fichier.close()

# Lecture du fichier U.txt
fichier = open("U.txt","r")
U = []
lignes = fichier.readlines()
for ligne in lignes :
	ligne = ligne.strip()
	ui = [float(elt) for elt in ligne.split("\n")]
	U.append(ui[0])
fichier.close()

# Lecture du fichier Uh.txt
fichier = open("Uh.txt","r")
Uh = []
lignes = fichier.readlines()
for ligne in lignes :
	ligne = ligne.strip()
	uhi = [float(elt) for elt in ligne.split("\n")]
	Uh.append(uhi[0])
fichier.close()

# Lecture du fichier Convergence.txt
fichier = open("Convergence.txt","r")
lignes = fichier.readlines()
NI, h, einf = [], [], []
for ligne in lignes :
	ligne = ligne.strip()
	ni, hi, einfi = [float(elt) for elt in ligne.split(",")]
	NI.append(int(ni))
	h.append(hi)
	einf.append(einfi)
fichier.close()


# ----------------------- GRAPHIQUES ----------------------- #
print("GRAPHIQUES\n")
# Graphiques pour le premier nombre d'intervalles
x = []
u = []
uh = []
for j in range(0,NI[0]+1): # On récupère les valeurs de x, u et uh pour le premier nombre d'intervalles
    x.append(X[j]) 
    u.append(U[j])
    uh.append(Uh[j])

err = [u_elt - uh_elt for u_elt, uh_elt in zip(u, uh)]

figure(1)
suptitle("Graphiques pour N = " + str(NI[0]) + "\n")
# Graphique des solutions
subplot(121)
plot(x,u,"b",label="solution exacte")
plot(x,uh,"r.", label="Solution approchee")
xlabel("Intervalle [O;L]")
legend()
title("Solutions")

# Graphique de l'erreur
subplot(122)
plot(x,err)
xlabel("Intervalle [O;L]")
title("Erreur U-Uh")
show()


# Graphiques pour le dernier nombre d'intervalles
x = []
u = []
uh = []
n = len(NI)
b = len(NI)
a = b - 1 
for i in range(0,a):
	a += NI[i]
b = a + NI[len(NI)-1] + 1

for j in range(a,b):
	x.append(X[j])
	u.append(U[j])
	uh.append(Uh[j])

err = [u_elt - uh_elt for u_elt, uh_elt in zip(u, uh)]


figure(2)
suptitle("Graphiques pour N = " + str(NI[len(NI)-1]) + "\n")
# Graphique des solutions
subplot(121)
plot(x,uh,"r.", label="Solution approchee")
plot(x,u,"b",label="solution exacte")
xlabel("Intervalle [O;L]")
legend()
title("Solutions")

# Graphique de l'erreur
subplot(122)
plot(x,err)
xlabel("Intervalle [O;L]")
title("Erreur U-Uh ")
show()


# ----------------------- CONVERGENCE ----------------------- #
from scipy import stats
x = -log(h)
y = -log(einf)
lr = stats.linregress(x, y)

figure(3)
plot(x, y, 'ro', label='erreurs sur les maillages')
plot(x, lr.intercept + lr.slope*x, 'b', label='moindres carrés')
x1 = x[len(x)-2]+0.5
y1 = lr.intercept + lr.slope*x1
x2 = x[len(x)-1]-0.5
y2 = lr.intercept + lr.slope*x2
plot([x1,x2,x2], [y1,y1, y2],'kD--', label="mesure de la pente")
legend()
p = lr.slope
text(x2+0.1,y1+0.1," p = {:.2f}".format(p),color='black');
xlabel("-log(h)")
ylabel("-log(einf)")
title("Courbe de convergence")
show()


print("FIN")

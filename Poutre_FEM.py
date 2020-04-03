#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Import des packages
import numpy as np
from numpy import linalg as lg
from matplotlib import pyplot as plt
 
# Parametres
E = 200e9 #Module de young
S = 0.01 # Section de la poutre
L = 1 # Longueur de la poutre
Fp = 1000 # Intensite de la force ponctuelle en Newtons
ne = 5 # Nombre delements
fl = 500 # Intensite de la force lineique en Newtons/metre
 
# Maillage
# Generation de la table coor
coor = np.zeros((ne+1,1)) # Initialisation
for i in range(ne+1): # Parcours des noeuds
    coor[i,0] = (float(i)/(ne))*L
#print (coor)
 
# Generation de la table connec
connec = np.zeros((ne,2))
for i in range(ne):
    connec[i,0] = i+1
    connec[i,1] = i+2
 
# Assemblage
K = np.zeros((ne+1,ne+1))
F = np.zeros((ne+1,1))
for i in range(ne):
    n1 = int(connec[i,0])
    n2 = int(connec[i,1])
    l = abs(coor[n1-1,0]-coor[n2-1,0])
    Ke = E*S/l*np.array([[1,-1],[-1,1]])
    K[n1-1,n1-1] += Ke[0,0]
    K[n2-1,n1-1] += Ke[1,0]
    K[n1-1,n2-1] += Ke[0,1]
    K[n2-1,n2-1] += Ke[1,1]
    F[n1-1] += 0.5*fl*l
    F[n2-1] += 0.5*fl*l

# Ajout de la force ponctuelle
F[ne] += Fp
 
# Prise en compte des CL
Kii = np.delete(K,0,0)
Kii = np.delete(Kii,0,1)
Fi = np.delete(F,0,0)
 
#plt.figure()
#plt.pcolor(Kii)
#plt.gca().invert_yaxis()
#plt.colorbar()
 
# Résolution
ui = lg.solve(Kii,Fi)
# On remet l'élément supprimé
u = np.insert(ui,0,0)
 
x = np.linspace(0,L,num=100)
uth = 1/(E*S)*((Fp+fl*L)*x-fl*np.square(x)/2)
 
plt.figure()
plt.plot(coor,u,x,uth)
plt.legend(('Deplacement element fini', 'Deplacement analytique'), loc=2)
plt.show()

print("Done.")
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
G2
"""
import numpy as np
from sympy import symbols,sqrt
import matplotlib.pyplot as plt

tamfig = (11,9)
tamlegend = 18
tamletra = 18
tamtitle = 18
tamnum = 16

a = symbols('a') #parametro de red
r = symbols('r') #distancia primeros vecinos

#%%

def energia_cohesion(epsilon,A_6,A_12):    
    num = epsilon*(A_6)**2
    den =   2*A_12  
    return -num/den

def energia_ion_par(EI,AE,energia_cohesion):
    rta = EI - AE - np.abs(energia_cohesion)
    return rta

def u_monio_Madelung(alpha,dist_segundos_vecinos,x):
    """
    Parameters
    ----------
    alpha : alpha (dato del problema)
    dist_segundos_vecinos : en funcion de la 
    distancia de 1 vecinos (que es r)
    x : r-/r+
    Raises
    ------
    TypeError
        error de dimension

    Returns
    -------
    u_monio = u sin e² y sin r_+
    para modelo de ion rigido (energia de Madelung)
    para cristales ionicos
    """
    cte = 2/dist_segundos_vecinos 
    xs = cte - 1
    rta1 = -alpha/(1+x)
    rta2 = -alpha/cte
    if x<= xs:
        rtaf = rta2
    else:
        rtaf = rta1
        
    return rtaf

def An(list_cant_vecinos, list_dist_vecinos,n):
    if len(list_cant_vecinos) != len(list_dist_vecinos):
        raise TypeError("Primera y segunda variable no tienen la misma long")
    rta = 0
    for k in range(len(list_cant_vecinos)):
        cant = list_cant_vecinos[k]
        dist = list_dist_vecinos[k]
        rta = rta + cant/(dist**n)
    return rta

def r_eq(sigma,A12,A6,list_cant_vecinos,list_dist_vecinos):
    A12 = An(list_cant_vecinos,list_dist_vecinos,12)
    A6 = An(list_cant_vecinos,list_dist_vecinos,6) 
    return sigma*((2*A12/A6)**(1/6))

def u_LJ(sigma,epsilon,list_cant_vecinos,list_dist_vecinos,r):
    A12 = An(list_cant_vecinos,list_dist_vecinos,12)
    A6 = An(list_cant_vecinos,list_dist_vecinos,6)
    return 2*epsilon*(A12*(sigma/r)**12 - A6*(sigma/r)**6)

#%% 

print('G2E1')


epsilon = 0.0031 #eV
sigma = 2.74 #Amnstrong

list_A6 = [8.40,12.2533,14.45392,14.45489]
list_A12 = [6.20, 9.11418 , 12.13188, 12.13229]
list_nombres = ['SC', 'BCC', 'FCC', 'HCP']

list_u = []
for j in range(4):
    A_6 = list_A6[j]
    A_12 = list_A12[j]
    u = energia_cohesion(epsilon,A_6,A_12)
    print(list_nombres[j], ':' , u)
    list_u.append(u)

mini = np.argmin(list_u)
print(list_nombres[mini])

#Estructura mas compacta ---> estructura mas estable

#%%

print('')
print('G2E4: energia de Madelung')

list_nombres = ['CsCl', 'NaCl', 'ZnCs']
list_colors = ['pink','blue','green']
list_dist_2_vecinos = [2/sqrt(3), sqrt(2), 4/sqrt(6)]
list_alphas = [1.7627, 1.7476, 1.6381]

list_x = np.linspace(0,1,500)
plt.figure(figsize=tamfig)
for j in range(3):
    label = list_nombres[j]
    dist_2_vecinos = list_dist_2_vecinos[j]
    alpha = list_alphas[j]
    list_u_monio = []
    for x in list_x:
        u = u_monio_Madelung(alpha,dist_2_vecinos,x)
        list_u_monio.append(u)
    plt.plot(list_x,list_u_monio,'-', lw = 3, color = list_colors[j],label = label)

plt.title('G2E4',fontsize = tamtitle)
plt.ylabel('u monio', fontsize=tamletra)
plt.xlabel('r-/r+', fontsize=tamletra)
plt.tick_params(labelsize = tamnum)
plt.legend(loc='best',markerscale=2,fontsize = tamlegend)

#%%

print('')
print('G2E5: panal de abejas')
# No es RB pero se puede usar las formulas porque phi solo depende de r (dist primeros vecinos)
# y en el panal de abejas r es cte

# 

sigma = symbols('sigma') 
epsilon = symbols('epsilon')

# interaccion hasta 2 vecinos
list_cant_vecinos = [3,6]
list_dist_vecinos = [a/a,a*sqrt(3)/a] #se normaliza por la dist a 1 vecinos, que es a

#u = u_LJ(sigma,epsilon,list_cant_vecinos,list_dist_vecinos,r)
#u_diff = Derivative(u, r)
#u_diff = u_diff.doit()
#r_eq = solve(Eq(u_diff,0),r)

A12 = An(list_cant_vecinos, list_dist_vecinos,12)
A6 = An(list_cant_vecinos, list_dist_vecinos,6)
r_eq_obt = r_eq(sigma,A12,A6,list_cant_vecinos,list_dist_vecinos)
u_eq = u_LJ(sigma,epsilon,list_cant_vecinos,list_dist_vecinos,r_eq_obt)

print('distancia de equilibrio entre 1 vecinos :', r_eq_obt)
print('energia de equilibrio:', u_eq)

#se puede evaluar:
sigma,epsilon = 2.25*1e-14,3.65

#%%

print('')
print('G2E6: modulo de bulk')


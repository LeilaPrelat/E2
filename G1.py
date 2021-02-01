# -*- coding: utf-8 -*-
"""
- Calcular los vectores de la RR (b) utilizando los vectores de la RD  (a)
escribiendolos en la base cartesiana

- G1E14

"""

import numpy as np
from sympy import symbols,sqrt,pi
from sympy import Symbol #vectores
from sympy import Matrix, Transpose #matrices
from sympy.physics.mechanics import dot,cross,ReferenceFrame

#%%
a = symbols('a') 
c = symbols('c') 

N = ReferenceFrame('N') #acceder a una base de vectores
x,y,z = Symbol('x'),Symbol('y'),Symbol('z') 

def Vector(coordenadas,dim):
    """
    Parameters
    ----------
    coordenadas : coordenadas de vector como array/list
    dim : dimension del vector (1,2 o 3)

    Raises
    ------
    TypeError
        error de dimension

    Returns
    -------
    vector : vector simbolico de dimension dim y 
            cuyas coord son la variable coordenadas

    """
    m = len(coordenadas)
    if m!=dim:
        raise TypeError(m,' = #coord != ', dim, '= dim del vector')
    if dim == 1:
        vector = coordenadas*N.x
    elif dim == 2 :
        cord1, cord2 = coordenadas
        vector = cord1*N.x + cord2*N.y 
    elif dim == 3:
        cord1,cord2,cord3 = coordenadas
        vector = cord1*N.x + cord2*N.y + cord3*N.z 
        
    return vector
    

def Volumen_CP(a1,a2,a3):
    """
    Parameters
    ----------
    a1 : 1 vector de la base VP
    a2 : 2 vector de la base VP
    a3 : 3 vector de la base VP
    Returns
    -------
    volumen de 1 CP
    """
    
    return dot(a1,cross(a2,a3))

def VP_RR(coord1,coord2,coord3):
    """
    Parameters
    ----------
    coord1 : coord del vector 1 VP
    coord2 : coord del vector 2 VP
    coord3 : coord del vector 3 VP
    Returns
    -------
    vectores primitivos de la red RR
    para 3D
    sirve para 2D: agregar + z en a1 y a2  
    para coord3 = [0,0,1], ignorar b3
    """
    
    dim1,dim2,dim3 = len(coord1),len(coord2),len(coord3)
    if dim1==dim2==dim3:
        dim = dim1
    else:
        raise TypeError('Error:Dimension de coordenadas diferentes')
    
    a1 = Vector(coord1,dim)
    a2 = Vector(coord2,dim)
    a3 = Vector(coord3,dim)
    vol = Volumen_CP(a1,a2,a3)
    
    b1 = 2*pi*cross(a2,a3)/vol
    b2 = 2*pi*cross(a3,a1)/vol
    b3 = 2*pi*cross(a1,a2)/vol
    return b1,b2,b3

def Matrix_coord(list_vectores):
    M = []
    for vector in list_vectores:
        coordx_v, coordy_v, coordz_v = dot(N.x,vector), dot(N.y,vector), dot(N.z,vector)
        M.append([coordx_v,coordy_v,coordz_v])
    return Matrix(M)

def vector_coord(vector):
    coordx_v, coordy_v, coordz_v = dot(N.x,vector), dot(N.y,vector), dot(N.z,vector)
    return [coordx_v, coordy_v, coordz_v]

def new_ind_miller(old_ind_miller,old_VPP_RD,new_VPP_RD,dim):
    """
    Parameters
    ----------
    old_ind_miller : [h,k,l] ind de miller de la base vieja
    old_VPP_RD : coordenadas de los viejos vectores primitivos de la red directa 
    (ai,aj,ak)
    new_VPP_RD : coordenadas de los nuevos vectores primitivos de la red directa 
    (a1,a2,a3), quiero los indices de miller en esta base
    dim : dimension de los VPP de la RD 
        
    Returns
    -------
    new_ind_miller : [h,k,l] ind de miller de los new_VPP_RD
    """
    coordi,coordj,coordk = old_VPP_RD
    # ai = Vector(coordi,dim)   
    # aj = Vector(coordj,dim)   
    # ak = Vector(coordk,dim)   

#   print('VP viejos de la RR:', bi,bj,bk)
    bi,bj,bk = VP_RR(coordi,coordj,coordk)
    
    coord1,coord2,coord3 = new_VPP_RD
    # a1 = Vector(coord1,dim)  
    # a2 = Vector(coord2,dim) 
    # a3 = Vector(coord3,dim) 

#   print('VP nuevos de la RR:', b1,b2,b3)    
    b1,b2,b3 = VP_RR(coord1,coord2,coord3)

    B_ijk = Matrix_coord([bi,bj,bk]) #base vieja
    B_xyz = Matrix_coord([b1,b2,b3]) #base nueva

#    print('Cambio de base')

    BT_xyz = Transpose(B_xyz) #transposer la matriz de coordenadas
    BT_ijk = Transpose(B_ijk)
    C_base = BT_xyz.inv()*BT_ijk
    C_base = np.matrix(C_base)

    ind_nuevo = np.array(C_base.dot(old_ind_miller))[0]

    return ind_nuevo

#%%

dim = 3
print('Ejemplo: hexagonal 3D (lado del hexaogo a y altura c)')

coord1 = np.array([a*0.5*sqrt(3),a*0.5,0])
a1 = Vector(coord1,dim)

coord2 = np.array([-a*0.5*sqrt(3),a*0.5,0])
a2 = Vector(coord2,dim)

coord3 = np.array([0,0,c])
a3 = Vector(coord3,dim)

#a1 = a*0.5*(sqrt(3)*N.x + 1*N.y + 0*N.z)
#a2 = a*0.5*(-sqrt(3)*N.x + 1*N.y + 0*N.z)
#a3 = c*(0*N.x + 0*N.y + 1*N.z)

vol = Volumen_CP(a1,a2,a3)
b1,b2,b3 =  VP_RR(coord1,coord2,coord3)

print('Volumen de 1CP:', vol)
print('VP de la RR:', b1,b2,b3)

#%%

print('')
print('Ejemplo: red SC de lado a')

coord1 = np.array([a,0,0])
coord2 = np.array([0,a,0])
coord3 = np.array([0,0,a])

a1 = Vector(coord1,dim)
a2 = Vector(coord2,dim)
a3 = Vector(coord3,dim)

vol = Volumen_CP(a1,a2,a3)
b1,b2,b3 = VP_RR(coord1,coord2,coord3)

print('Volumen de 1CP:', vol)
print('VP de la RR:', b1,b2,b3)
print('La RR de una SC de lado a es una SC de lado 2*pi/a')

#%%

print('')
print('Ejemplo: red FCC de lado a')

coord1 = np.array([a*0.5,0,a*0.5])
coord2 = np.array([a*0.5,a*0.5,0])
coord3 = np.array([0,a*0.5,a*0.5])

a1 = Vector(coord1,dim)  
a2 = Vector(coord2,dim)   
a3 = Vector(coord3,dim)  

vol = Volumen_CP(a1,a2,a3)
b1,b2,b3 =  VP_RR(coord1,coord2,coord3)

print('Volumen de 1CP:', vol)
print('VP de la RR:', b1,b2,b3)
print('La RR de una FCC de lado a es una BCC de lado 2*pi/a')

#%%

print('')
print('Ejemplo: red BCC de lado a')

coord1 = np.array([a*0.5,0,0])
coord2 = np.array([0,a*0.5,0])
coord3 = np.array([a*0.5,a*0.5,a*0.5])

a1 = Vector(coord1,dim)   
a2 = Vector(coord2,dim) 
a3 = Vector(coord3,dim) 

vol = Volumen_CP(a1,a2,a3)
b1,b2,b3 =  VP_RR(coord1,coord2,coord3)

print('Volumen de 1CP:', vol)
print('VP de la RR:', b1,b2,b3)
print('La RR de una BCC de lado a es una FCC de lado 2*pi/a')

#%%

print('')
print('G1E14: red FCC de lado a')

coordi = np.array([0,a*0.5,a*0.5])
coordj = np.array([a*0.5,0,a*0.5])
coordk = np.array([a*0.5,a*0.5,0])

old_VPP_RD = [coordi, coordj, coordk]

coord1 = np.array([a,0,0])
coord2 = np.array([a*0.5,a*0.5,0])
coord3 = np.array([0,a*0.5,a*0.5])

new_VPP_RD = [coord1,coord2,coord3]

ind1 = [1,0,0] #h, k l
ind2 = [0,0,1] #h, k l

ind1_nuevo = new_ind_miller(ind1,old_VPP_RD,new_VPP_RD,dim)
ind2_nuevo = new_ind_miller(ind2,old_VPP_RD,new_VPP_RD,dim)
print(ind1_nuevo,ind2_nuevo)

#%%

#K1 = np.dot(ind1,[bi,bj,bk]) # K = h*bi + k*bj + l*bk
#K2 = np.dot(ind2,[bi,bj,bk])

#K1 == -b1 + b3   # --- > (-1,0,1)
#K2 == b1 + b2     # --- > (1,1,0)

#%%

print('')
print('G1E15')

def factor_estructura(h,k,l):
    return 1 + (-1)**(h+k) + (-1)**(h+l) + (-1)**(l+k)

Sk = []
for h in [0,1]:
    for k in [0,1]:
        for l in[0,1]:
            Sk.append(factor_estructura(h,k,l))

print(set(Sk))


#%%





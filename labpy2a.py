# -*- coding: utf-8 -*-
"""
Spyder Editor
"""

#Boundary value problem for  Poisson Equation with Dirichlet boundary condition
#Reducing to a SLAE
import numpy as np
import math as M
import matplotlib.pyplot as plt
from scipy import linalg
#******************************************************************************

#Initialisation of variables
a = 4
b = 4
n = 4
m = 4
h1 = a/(n + 1)
h2 = b/(m + 1)
alpha = (h1/h2)**2

x = np.zeros(n + 2)       
y = np.zeros(m + 2)
#f = np.zeros((n + 2, m + 2))
u = np.zeros((n + 2, m + 2))
uex = np.zeros((n + 2, m + 2))
    
#Setting Mesh
for i in range (0, n + 2):
        x[i] = i*h1
for j in range (0, m + 2):
        y[j] = j*h2
    
#Setting Boundary conditions for y:
for j in range (0, m + 2):
        u[n + 1, j] = M.sinh(y[j])/M.sinh(b) 
        
#Setting Boundary conditions for x:
for i in range (0, n + 2):
        u[i, m + 1] = M.sin(x[i])/M.sin(a)   
       
#for i in range (0,n,1):
#    z[i] = i
#xa = np.arange(n+2)*h1
    
#Numerical Solution (SLAE method):
nm=n*m

Ad = np.zeros((m,m))
Au = np.zeros((m,m))
Az = np.zeros((m,m))
A = np.zeros((nm,nm))


i,j=np.indices(Ad.shape)
i,j=np.indices(Au.shape)
Ad[i==j]=-2.*(1.+alpha)
Ad[i==j+1]=alpha
Ad[i+1==j]=alpha
Au[i==j] = 1.



for i in range(n):
    for i1 in range(m):
        for i2 in range(m):
             A[i1+i*m,i2+i*m]=Ad[i1,i2]
 
for i in range(n-1):
   for i1 in range(m):
        for i2 in range(m):
             A[i1+(i+1)*m,i2+i*m]=Au[i1,i2]
             
for i in range(n-1):
   for i1 in range(m):
        for i2 in range(m):
             A[i1+i*m,i2+(i+1)*m]=Au[i1,i2]
             
#Inhomogeneous term  

print(A)
ft = np.zeros((n, m))

for i in range (n):
    for j in range (m):
            ft[i, j] = x[i+1]*y[j+1]*0    

for j in range (m):
    ft[0, j] = ft[0, j] - u[0, j+1]
    ft[n-1, j] = ft[n-1, j] - u[n + 1, j+1]
        
for i in range (n):
    ft[i, 1] = ft[i, 1] - alpha*u[i+1, 0]
    ft[i, m-1] = ft[i, m-1] - alpha*u[i+1, m + 1]
    
F = np.zeros(nm)
for i in range (n):
    for j in range (m):
        k = i*m + j
        F[k] = ft[i, j]

U = np.zeros(nm)
F = F.transpose()

U = linalg.solve(A,F)
#u = U.reshape(n,m)
    #for i in range (nm):
        #print (F[i])

for i in range (n):
    for j in range (m): 
            #k = i*m + j
        u[i + 1, j + 1] = U[i*m + j]

#Exact Solution:
for i in range (0, n + 2):
    for j in range (0, m + 2):
            uex[i, j] = M.sin(x[i])*M.sinh(y[j])/(M.sin(a)*M.sinh(b))
            
fig, ax = plt.subplots()

plt.plot(x,u[:,3]) ### -2: second last element
plt.plot(x,uex[:,3],ls='dashed')
#plt.plot(ya,u[-2,:])
#plt.plot(ya,ulin[-2,:],ls='dashed')
plt.show()
fig, ax = plt.subplots()

plt.plot(y,u[3,:]) ### -2: second last element
plt.plot(y,uex[3,:],ls='dashed')
#plt.plot(ya,u[-2,:])
#plt.plot(ya,ulin[-2,:],ls='dashed')
plt.show()

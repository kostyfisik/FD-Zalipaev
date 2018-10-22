# -*- coding: utf-8 -*-
"""
Spyder Editor
"""

#This program solves Poisson Equation with Iteration Method

import numpy as np
import math as M
import matplotlib.pyplot as plt
#******************************************************************************

#Initialisation of variables
a = 2
b = 2
n = 20
n1=n/2
print(n)
m = 20
m1=m/2
h1 = a/(n + 1)
h2 = b/(m + 1)
H = 1/(2/(h1*h1) + 2/(h2*h2))
N = 200

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
        
#Setting Boundarry conditions  for x:
for i in range (0, n + 2):
        u[i, m + 1] = M.sin(x[i])/M.sin(a)      
    
#Numerical Solution (iteration method):
for k in range (N):
    for i in range (1, n+1):
        for j in range (1, m+1):
                u[i, j] = H*((u[i - 1, j]
                + u[i + 1, j])/(h1*h1) + 
                 (u[i, j - 1] + u[i, j + 1])/(h2*h2))
#print (u)

#Exact Solution:
for i in range (0, n + 2):
    for j in range (0, m + 2):
            uex[i, j] = M.sin(x[i])*M.sinh(y[j])/(M.sin(a)*M.sinh(b))
            
# Plotting solution            
fig, ax = plt.subplots()

plt.plot(x,u[:,10]) ### -2: second last element
plt.plot(x,uex[:,10],ls='dashed')
#plt.plot(ya,u[-2,:])
#plt.plot(ya,ulin[-2,:],ls='dashed')
plt.show()

fig, ax = plt.subplots()

plt.plot(y,u[10,:]) ### -2: second last element
plt.plot(y,uex[10,:],ls='dashed')
#plt.plot(ya,u[-2,:])
#plt.plot(ya,ulin[-2,:],ls='dashed')
plt.show()
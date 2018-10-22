# -*- coding: utf-8 -*-
"""
Created on Sat May  5 19:32:25 2018

@author: victor
"""
#Wave equation boundary value problem: Explicit method
import numpy as np
import math as M
import matplotlib.pyplot as plt
#*********************************************************************************

#Initialisation of variables
pi=3.141592653589793
T = 1.5
n = 80
m = 176
h = 1/(n+1)
k = T/m
s = (k/h)**2
A=np.zeros((n,n))
u=np.zeros((m+1,n+2))
uex=np.zeros((m+1,n+2))
U=np.zeros(n)
V=np.zeros(n)
W=np.zeros(n)
error=np.zeros(m+1)

#Setting Mesh
x = h*np.arange(n+2)
t = k*np.arange(m+1)
#Exact Solution:
for i in range (0, n + 2):
    for j in range (0, m + 1):
            uex[j, i] = M.sin(pi*x[i])*(M.cos(pi*t[j])+M.sin(pi*t[j])/pi)
     
x[0]=0.0
x[n+1]=1.0
for i in range (0, n):
    W[i] = M.sin(pi*h*(i+1))
    u[0,i+1] = W[i]
    V[i] = M.sin(pi*h*(i+1))*k+M.sin(pi*h*(i+1))
    u[1,i+1] = V[i]
for j in range (0, m + 1):
    u[j,0]=0.0
    u[j,n+1]=0.0
    
#Explicit method
i,j=np.indices(A.shape)
A[i==j]=2-2.0*s
A[i==j+1]=s
A[i+1==j]=s

for j in range (2, m + 1):
    U=A.dot(V)-W
    W=V
    V=U
    for i in range (0, n):
        u[j,i+1] = U[i]
    
plt.plot(x,u[m,:]) ### -2: second last element
plt.plot(x,uex[m,:],ls='dashed')
plt.show()

for j in range (0, m + 1):
    sum=0.0
    for i in range (0, n + 2):
        sum=sum+(uex[j,i]-u[j,i])**2
    error[j]=M.sqrt(sum/(n+2))
plt.plot(t,error) ### -2: second last element
plt.show() 

fig1 = plt.figure(1)
ColorMap(uex, 'Exact Solution_Good')
plt.savefig("Exact.png",pad_inches=0.02, bbox_inches='tight')

plt.show()

fig2 = plt.figure(2)
ColorMap(u, 'Numerical Solution_Good')
plt.savefig("Numerical.png",pad_inches=0.02, bbox_inches='tight')

plt.show() 
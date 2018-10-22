# -*- coding: utf-8 -*-
"""
Created on Sat May  5 19:32:25 2018

@author: victor
"""
#Heat transfer boundary value problem: Explicit, implicit, Crank-Nocolson methods
import numpy as np
import math as M
import matplotlib.pyplot as plt
#*********************************************************************************

#Initialisation of variables
pi=3.141592653589793
T = 0.5
n = 20
m = 500
h = 1/(n+1)
k = T/m
s = k/h/h
A=np.zeros((n,n))
Id=np.zeros((n,n))
B=np.zeros((n,n))
u=np.zeros((m+1,n+2))
unum=np.zeros((m+1,n+2))
uex=np.zeros((m+1,n+2))
V0=np.zeros(n)
V=np.zeros(n)
error=np.zeros(m+1)

#Setting Mesh
x = h*np.arange(n+2)
t = k*np.arange(m+1)
#Exact Solution:
for i in range (0, n + 2):
    for j in range (0, m + 1):
            uex[j, i] = M.sin(pi*x[i])*M.exp(-pi**2*t[j])
     
i,j=np.indices(B.shape)
i,j=np.indices(B.shape)
B[i==j]=2.0
B[i==j+1]=-1.0
B[i+1==j]=-1.0
Id[i==j] = 1.0

for i in range (0, n):
    V0[i] = M.sin(pi*h*(i+1))
    u[0,i+1] = V0[i]
for j in range (0, m + 1):
    u[j,0]=0.0
    u[j,n+1]=0.0
    
#Explicit method
i,j=np.indices(A.shape)
A[i==j]=1-2.0*s
A[i==j+1]=s
A[i+1==j]=s

#V=V0.transpose()
V=np.transpose(V0)

#for i in range (n):
#    print (V[i])
    
#print (V)

for j in range (1, m + 1):
#    V=A.dot(V)
    V=np.dot(A,V)
    for i in range (0, n):
        u[j,i+1] = V[i]
        
#unum=u
    
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

#%Implicit method

i,j=np.indices(A.shape)
A[i==j]=1+2.0*s
A[i==j+1]=-s
A[i+1==j]=-s

V=V0.transpose()

for j in range (1, m + 1):
    V=Lin.inv(A).dot(V)
    for i in range (0, n):
        u[j,i+1] = V[i]
        
#unum=u
    
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

#Crank-Nocolson method

V=V0.transpose()


A=Lin.inv(2*Id+s*B).dot(2*Id-s*B)
for j in range (1, m + 1):
    V=A.dot(V)
    for i in range (0, n):
        u[j,i+1] = V[i]
        
unum=u
    
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
ColorMap(unum, 'Numerical Solution_Good')
plt.savefig("Numerical.png",pad_inches=0.02, bbox_inches='tight')

plt.show()
 

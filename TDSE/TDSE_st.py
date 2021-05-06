#!/usr/bin/python3.7

import numpy as np
import matplotlib.pyplot as plt
L = 30.0
x0 = -5.0
sig = 0.5
dx = 0.5
dt = 0.02
k = 1.0
w=2
K=w**2
a=np.power(K,0.25)
xs = np.arange(-L,L,dx)
nn = len(xs)

mu = k*dt/(dx)**2
dd = 1.0+mu
ee = 1.0-mu
ti = 0.0
tf = 100.0
t = ti
V=np.zeros(len(xs))
u=np.zeros(nn,dtype="complex") 
V=K*(xs)**2/2            #harmonic oscillator potential
u=(np.sqrt(a)/1.33)*np.exp(-(a*(xs - x0))**2)+0j    #initial condition for wave function
u[0]=0.0          #boundary condition
u[-1] = 0.0      #boundary condition

A = np.zeros((nn-2,nn-2),dtype="complex")     #define A
for i in range(nn-3):
    A[i,i] = 1+1j*(mu/2+w*dt*xs[i]**2/4)
    A[i,i+1] = -1j*mu/4.
    A[i+1,i] = -1j*mu/4.
A[nn-3,nn-3] = 1+1j*mu/2+1j*dt*xs[nn-3]**2/4 

B = np.zeros((nn-2,nn-2),dtype="complex")    #define A*
for i in range(nn-3):
    B[i,i] = 1-1j*mu/2-1j*w*dt*xs[i]**2/4
    B[i,i+1] = 1j*mu/4.
    B[i+1,i] = 1j*mu/4.
    B[nn-3,nn-3] = 1-1j*(mu/2)-1j*dt*xs[nn-3]**2/4

X = np.linalg.inv(A)    #take inverse of A
plt.ion()
l, = plt.plot(xs,np.abs(u),lw=2,color='blue')   #plot initial wave function
T=np.matmul(X,B)                                #multiply A inverse with A*

while t<tf:
    u[1:-1]=np.matmul(T,u[1:-1]) #updating u but leaving the boundary conditions unchanged
    l.set_ydata((abs(u)))              #update plot with new u
    t += dt
    plt.pause(0.00001)

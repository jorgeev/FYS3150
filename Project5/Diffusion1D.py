#! /usr/bin/python
# -*- coding: utf-8 -*-
# Diffusion1D.py

# Version: 2018.12.14

# Requiered libraries
import sys, os, math, time
import numpy as np
import scipy.linalg as cp
import matplotlib.pyplot as plt
from numba import jit, njit, prange

@njit() 
def v(x,t,n): # Analytic solution
    Bn = (n*math.pi*math.cos(n*math.pi) - math.sin(n*math.pi))/((n*math.pi)**2)
    serie = math.sin(n*math.pi*x)*math.exp(-t*(n*math.pi)**2)
    return Bn*serie

def cn(alpha,u,N,T): #Crank-Nicolson Method
    N = N+1
    aux = np.zeros(N)
    aux[0] = 1.+alpha
    aux[1] = -alpha/2
    A = cp.toeplitz(aux)
    A[0] = A[-1]= 0 # Set first and last row to 0
    A[0,0] = A[-1,-1] = 1
    A = cp.inv(A, overwrite_a=True, check_finite=True)
    for tt in range(1,T):
        forward(alpha/2,u[tt],u[tt-1],N) # Explicit scheme
        #explicit(N,T,u[tt],alpha/2)
        u[tt] = A@u[tt]
    return(u)

def forward(alpha,u,uPrev,N): # Explicit scheme
    betta = 1.-2.*alpha
    for x in range(1,N-1): #loop from i=1 to i=N
        u[x] = alpha*uPrev[x-1] + betta*uPrev[x] + alpha*uPrev[x+1]
        
def implicit(alpha,u,N,T): #Implicit method
    N = N+1
    aux = np.zeros(N)
    aux[0] = 1.+2.*alpha
    aux[1] = -alpha
    A = cp.toeplitz(aux)
    A[0] = A[-1]= 0 # Set first and last row to 0
    A[0,0] = A[-1,-1] = 1 # set values to considere boundary conditons different than 0
    A = cp.inv(A, overwrite_a=True, check_finite=True)
    for ii in range(1,T):
        uaux = u[ii-1]
        tt = A@uaux
        u[ii] = tt
    return(u)

@njit()
def explicit(X,T,u,alpha):
    b = 1.-2.*alpha
    for ii in range(0,T-1):
        for jj in range(1,X):
            u[ii+1,jj] = alpha*u[ii,jj-1] + (b)*u[ii,jj] + alpha*u[ii,jj+1]
    return u     

# Iniciate conditions

N = 100 
dt = 5e-5
T = int(0.25/dt)+1 # Target final time
mode = 2 # alterneate 0 = explicit, 1 = implicit, 2 = CN

u = np.zeros((T,N+1), np.double)
(x,dx) = np.linspace(0,1,N+1, retstep=True)

# Iniciate vector to store analytic and method solution
us = np.zeros((2,N+1), np.double)
ua = np.zeros((2,N+1), np.double)
nx = 100 # Sum series counter

for xx in range(N+1): # Calcualate analytic solution at two different times
    for nn in range(1,nx+1):
        ua[0,xx]+=v(x[xx],0.15,nn)
    ua[0,xx] *= 2
    ua[0,xx] += x[xx]
    
for xx in range(N+1):
    for nn in range(1,nx+1):
        ua[1,xx]+=v(x[xx],0.25,nn)
    ua[1,xx] *= 2
    ua[1,xx] += x[xx]
    
alpha = dt/dx**2

u[0,0] = 0. #Boundary conditions
u[:,-1]= 1.


# Calculations
if mode == 0:
    t0 = time.time()
    u = explicit(N,T,u,alpha)
    t1 = time.time()
    print(t1-t0)
elif mode == 1:
    t0 = time.time()
    u = implicit(alpha,u,N,T)
    t1 = time.time()
    print(t1-t0)
elif mode == 2:
    t0 = time.time()
    u = cn(alpha,u,N,T)
    t1 = time.time()
    print(t1-t0)
else:
    print("ERROR")




y = dt*np.linspace(0,T,T+1) # Vector y for ploting
plt.plot(x,ua[0], '-', label = "T= 0.15", color='red') #analytic solutions
plt.plot(x,ua[1], '--', label = "T= 0.25", color='red')

# Look for speciffic solutions plot them and store result for error calculations
for ii in range(0,T+1):
    if y[ii] == 0.15:
        us[0,:] = u[ii,:]
        plt.plot(x,u[ii,:], ':', label = "T= "+str(y[ii]), color='C0')
    elif y[ii] == 0.25:
        us[1,:] = u[ii,:]
        plt.plot(x,u[ii,:], ':', label = "T= "+str(y[ii]), color='C1')
plt.legend(loc='best')
plt.xlabel('Length [L]')
plt.ylabel('Temperature')
plt.tight_layout()
#plt.show()
#plt.savefig("Plot.pdf", format='pdf')
del u


# Compute relative error and save files
err = np.zeros((2,N+1), np.double)
err[:,1:] = (ua[:,1:]-us[:,1:])/ua[:,1:]

f1 = 'Method15'
f2 = 'Method25'
x = x.reshape(x.shape+(1,))
ua15 = ua[0].reshape(ua[0].shape+(1,))
us15 = us[0].reshape(us[0].shape+(1,))
err0 = err[0].reshape(err[0].shape+(1,))
Aux = np.hstack((x,ua15,us15, err0))
np.savetxt(f1, Aux, fmt='%0.6e', delimiter=',')
ua15 = ua[1].reshape(ua[1].shape+(1,))
us15 = us[1].reshape(us[1].shape+(1,))
err0 = err[1].reshape(err[1].shape+(1,))
Aux = np.hstack((x,ua15,us15, err0))
np.savetxt(f2, Aux, fmt='%0.6e', delimiter=',')


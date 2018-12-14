import math, sys, os
import numpy as np
import scipy.linalg as cp
from  matplotlib import pyplot as plt
from numba import njit

@njit
def jacobiHE2D(U, alpha, N, T):
    r = 1./(1.+4.*alpha)
    Atemp = U.copy() #Reference matrix
    for xx in range(1,N): #Iterate over the whole mesh
        for yy in range(1,N):
            U[xx,yy] = r*(alpha*(Atemp[xx,yy+1]+Atemp[xx,yy-1]+Atemp[xx+1,yy]+Atemp[xx-1,yy])+Atemp[xx,yy])
    return U



N = 100 # Square mesh size
x,dx = np.linspace(0,1,N+1,retstep=True)
y = x.copy() #Square mesh

dt = 0.01
Tmax = 1.
T = int(Tmax/dt)+1 #Target max time

alpha = dt/dx**2

#Initialize Mesh
U = np.zeros((N+1,N+1))

#Saving initial plot
fig, ax = plt.subplots()
im = ax.pcolormesh(x,y,U,cmap='jet',vmin=0., vmax=1.)
cbar = fig.colorbar(im)
plt.savefig("test0"+".png", format='png')

# Set Boundary conditions
U[:,0]  = 0.
U[:,-1] = 1.
U[0,:]  = 1.
U[-1,:] = 0.

print("Number of iterations to reach T_max =",Tmax,": ", T)

img = 1
for tt in range(1,T):
	U = jacobiHE2D(U, alpha, N, T)
	# Start plotting
	fig, ax = plt.subplots()
	im = ax.pcolormesh(x,y,U,cmap='jet',vmin=0., vmax=1.)
	cbar = fig.colorbar(im)
	plt.title("T ="+str(dt*tt))
	plt.savefig("test"+str(img)+".png", format='png')
	del fig,ax,im,cbar # Liberate Memory
	img += 1


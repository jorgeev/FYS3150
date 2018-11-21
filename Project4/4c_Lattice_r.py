#! /usr/bin/python
# -*- coding: utf-8 -*-
# 4cLattice20.py

# Version: 2018.11.11.01


# Requiered libraries
import numpy as np
from numpy.random import random as nrand
import scipy.linalg as cp
import time
import random
import math
import matplotlib.pyplot as plt
from numba import jit
import time

# Monte Carlo Alg
@jit(nopython=True)
def MC(temperature, spins, MCc, MSpins): # Temperature, Number of Spins, Number of cycles
    PBC = lambda idx,lim,add: (idx + lim + add) % lim

    
    E = 0.; Eavr = 0.
    M = 0.; Mavr = 0.
    lMl = 0.
    VE = np.zeros(MCc+1)
    VM = np.zeros(MCc+1)
    
    # Sampling rule
    w = np.zeros(17,np.float64)
    for de in range(-8,9,4):
        w[de+8] = math.exp(-de/temperature)
    
    # Initial energy and magnetization
    M = MSpins.sum()
    for ii in range(spins):
        for jj in range(spins):
            E -= MSpins[ii,jj] * (MSpins[PBC(ii,spins,-1),jj] + MSpins[ii,PBC(jj,spins,1)])
    
    VE[0] = E
    VM[0] = M
    
    #MC
    for ii in range(MCc):
        # print(ii)
        for ss in range(spins**2):
            x = int(nrand() * spins)
            y = int(nrand() * spins)
            # print(x,y)
            dE = 2 * MSpins[x,y] * (MSpins[PBC(x,spins,-1),y] + MSpins[PBC(x,spins,1),y] + MSpins[x,PBC(y,spins,-1)] + MSpins[x, PBC(y,spins,1)])
            # print(dE)
            if nrand() <= w[dE+8]:
                MSpins[x,y] *= -1
                M += 2 * MSpins[x,y]
                E += dE
                
              
        Eavr += E;
        VE[ii+1] = Eavr/(ii+1)
        lMl += int(math.fabs(M))
        VM[ii+1] = lMl/(ii+1)

    #Eavr /= float(MCc*  spins**2)
    #Mavr /= float(MCc)
    #lMl /= float(MCc * spins**2)
    VE /= float(spins**2)
    VM /= float(spins**2)

    return (VE, VM)


t0 = 1.; ti = 2.4; steps = 2
nspins = 20; cycls = 1e6;
fname = "random"

Temps = np.linspace(t0,ti,steps)

#MSpins = np.ones((nspins,nspins),np.int8) # Spins Matrix pointing up
MSpins = np.random.rand(nspins, nspins)
for ii in range(nspins):
    for jj in range(nspins):
        if MSpins[ii,jj] >= 0.5:
            MSpins[ii,jj] = 1
        else:
            MSpins[ii,jj] = -1
MSpins = np.int_(MSpins)
#print(MSpins)

Energy = np.zeros([steps,int(cycls)+1]); Magnetization = np.zeros([steps,int(cycls)+1])

for stp in range(steps):
    Energy[stp], Magnetization[stp]= MC(Temps[stp], nspins, int(cycls), MSpins)


Aux = np.hstack((np.transpose(Energy), np.transpose(Magnetization)))
np.savetxt(fname, Aux, fmt='%0.6f', delimiter='   ')


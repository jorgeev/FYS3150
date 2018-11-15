#! /usr/bin/python
# -*- coding: utf-8 -*-
# Dirichlet.py

# Version: 2018.11.4.01


# Requiered libraries
from __future__ import division
import time, random, math, time, sys, os
import matplotlib.pyplot as plt
from numba import jit, prange, njit, int32
import numpy as np
from numpy.random import random as nrand

@njit(int32(int32, int32, int32)) # Periodic Boundary Condition
def PBC(idx, lim, add):
    return (idx + lim + add) % lim

# Monte Carlo Alg
def MC(temperature, spins, MCc): # Temperature, Number of Spins, Number of cycles
    MSpins = np.ones((spins,spins), np.int32) # Spins Matrix pointing up
    
    # Initialize Statistics
    E = 0.
    Eavr = 0.; Evar = 0.; E2av = 0.
    Mavr = 0.; Mvar = 0.; M2av = 0.
    lMl = 0.

    # Initial Magnetization and Energy
    M = MSpins.sum()
    for ii in range(spins):
        for jj in range(spins):
            E -= MSpins[ii,jj] * (MSpins[PBC(ii,spins,-1),jj] + MSpins[ii,PBC(jj,spins,1)])
    
    # Sampling rule
    w = np.zeros(17,np.float64)
    for de in range(-8,9,4):
        w[de+8] = math.exp(-de/temperature)
    
    #Monte Carlo Computation
    Eavr, E2av, Mavr, M2av, lMl = nMCC(MCc,spins, MSpins, E, M, Eavr, E2av, Mavr, M2av, lMl, w)

    Eavr /= float(MCc); E2av /= float(MCc)
    Mavr /= float(MCc); M2av /= float(MCc)
    lMl /= float(MCc * spins**2)
    
    # Calculate Variance
    E_variance = (E2av - (Eavr**2)) / float(spins**2)
    #M_variance = (M2av - (Mavr**2)) / float(spins**2)
    M_variance = (M2av - ((lMl* spins**2)**2)) / float(spins**2)
    SH = E_variance / float(temperature**2)
    SS = M_variance / float(temperature)

    Eavr /= float(spins**2)
    Mavr /= float(spins**2)
    E2av /= float(spins**2)
    M2av /= float(spins**2)

    #return (Eavr, E_variance, Mavr, M_variance, lMl)
    return (Eavr, SH, Mavr, SS, lMl, E_variance, M_variance)

def MCC(MCc, spins, MSpins, E, M, Eavr, E2av, Mavr, M2av, lMl, w):
    for ii in range(MCc):
        for ss in range(spins**2):
            x = int(nrand() * spins) #Pick up a spin
            y = int(nrand() * spins) #Pick up a spin
            dE = 2 * MSpins[x,y] * (MSpins[PBC(x,spins,-1),y] + MSpins[PBC(x,spins,1),y] + MSpins[x,PBC(y,spins,-1)] + MSpins[x, PBC(y,spins,1)])
            
            # Metropolis
            if nrand() <= w[dE+8]:
                MSpins[x,y] *= -1
                M += 2 * MSpins[x,y]
                E += dE
                
        Eavr += E; E2av += E**2
        Mavr += M; M2av += M**2
        lMl += int(math.fabs(M))
    return(Eavr, E2av, Mavr, M2av, lMl)


def launcher(Energy, S_Heat, Magnetization, Susceptibility, lMagnetizationl, Ev, Mv, Temps, nspins, cycls, steps):
    for stp in prange(steps):
        Energy[stp],S_Heat[stp],Magnetization[stp],Susceptibility[stp],lMagnetizationl[stp], Ev[stp], Mv[stp] = nMC(Temps[stp], nspins, int(cycls))
    return(Energy, S_Heat, Magnetization, Susceptibility, lMagnetizationl, Ev, Mv)

# Run parameters
if __name__== "__main__":
    args = sys.argv
    if len(args) > 5:
        print("Using command line inputs:")
        name = "Lattice_" + args[1]
        nspins = int(args[1])
        t0 = float(args[2])
        ti = float(args[3])
        steps = int(args[4])+1
        cycls = float(args[5])       
    else:
        print("Demo mode")
        t0 = 2.; ti = 3.; steps = 10+1
        nspins = 2; cycls = 1e6;
        name = "Demo"


sline = str('Lattice: {0}^2, T[{1:0.1f},{2: 0.1f}], #Steps = {3:0}, #Cycles = {4:0.0f} '.format(nspins, t0, ti, steps, cycls))

print(sline)

# Vector with all temps to be computed
Temps = np.linspace(t0,ti,steps)

#Precompile all functions
nl = jit(launcher, nopython = True, parallel = True)
nMC = njit(MC)
nMCC = njit(MCC)

# Vectros for results
Energy = np.zeros(steps); Magnetization = np.zeros(steps)
Ev = np.zeros(steps); Mv = np.zeros(steps)
S_Heat = np.zeros(steps); Susceptibility = np.zeros(steps)
lMagnetizationl = np.zeros(steps)


c0 = time.time()

#Energy, S_Heat, Magnetization, Susceptibility, lMagnetizationl, Ev, Mv = nl(Energy, S_Heat, Magnetization, Susceptibility, lMagnetizationl, Ev, Mv, [Temps[0]], nspins, 1, 1)
#print("Done")

Energy, S_Heat, Magnetization, Susceptibility, lMagnetizationl, Ev, Mv = nl(Energy, S_Heat, Magnetization, Susceptibility, lMagnetizationl, Ev, Mv, Temps, nspins, cycls, steps)

c1 = time.time()

print("Time elapsed: ", c1-c0, "sec.\n")


# Temp       Energy      Variance   Sp_Heat    Magneti     absMagne   Variance  Suseptiblity
aux = np.hstack((Temps.reshape(Temps.shape+(1,)),Energy.reshape(Energy.shape+(1,)), Ev.reshape(Ev.shape+(1,)),S_Heat.reshape(S_Heat.shape+(1,))))
aux = np.hstack((aux, Magnetization.reshape(Magnetization.shape+(1,)),lMagnetizationl.reshape(lMagnetizationl.shape+(1,)),Mv.reshape(Mv.shape+(1,)),Susceptibility.reshape(Susceptibility.shape+(1,))))
np.savetxt(name, aux, fmt='%0.6f', delimiter='   ')


# And finally plot
f = plt.figure(figsize=(18, 10)); # plot the calculated values    

sp =  f.add_subplot(2, 2, 1 );
plt.plot(Temps, Energy, 'b-', color="green");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Energy ", fontsize=20);

sp =  f.add_subplot(2, 2, 2 );
plt.plot(Temps, abs(lMagnetizationl), 'b-', color="red");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Magnetization ", fontsize=20);

sp =  f.add_subplot(2, 2, 3 );
plt.plot(Temps, S_Heat, 'b-', color="blue");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Specific Heat ", fontsize=20);

sp =  f.add_subplot(2, 2, 4 );
plt.plot(Temps, Susceptibility, 'b-', color="black");
plt.xlabel("Temperature (T)", fontsize=20);
plt.ylabel("Susceptibility", fontsize=20);

plt.savefig(name+".png")
plt.show()




#plt.plot(Ex)
# np.savetxt("Ex.txt", Ex, fmt='%0.6f')
exit()

# Quantum dots
# Libraries
import numpy as np
import math
import scipy.linalg as cpy
from jtools import *
# import time

def dots(size = 10, RMin = 0.,RMax = 10.,lOrbital = 0):
    # Creates a tridiagonal matrix
    step = RMax / (size + 1.)
    Aux = np.zeros(size)
    Aux[0] = 2. / step**2
    Aux[1] = -1. / step**2
    toe = cpy.toeplitz(Aux)
    
    # Adding potential values
    OFactor = lOrbital * (lOrbital + 1)
    r = np.linspace(RMin,RMax, size)
    v = np.zeros(size)
    for ii in range(size):
        r[ii] = RMin + (ii+1) * step
        v[ii] = r[ii]**2 + OFactor / r[ii]**2
        toe[ii,ii] = toe[ii,ii] + v[ii]
    del Aux, r, v
    return(toe, size)

A, m = dots(10)
lam, V = jacobiT(A,m)
del A
print(lam)


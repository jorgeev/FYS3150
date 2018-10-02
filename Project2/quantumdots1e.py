# Quantum dots
# Libraries
import numpy as np
import math
import scipy.linalg as cpy
from jtools2 import *
# import time

r_error = lambda x, y : abs( (x - y) / x)

def dots(size = 5, RMin = 0.,RMax = 10.):
    step = RMax / (size + 1.)
    Aux = np.zeros(size)
    Aux[0] = 2. / step**2
    Aux[1] = -1. / step**2
    toe = cpy.toeplitz(Aux)
    
    # Adding potential values
    r = np.linspace(RMin,RMax, size)
    v = np.zeros(size)
    for ii in range(size):
        r[ii] = RMin + (ii+1) * step
        v[ii] = r[ii]**2
        toe[ii,ii] = toe[ii,ii] + v[ii]
    del Aux, r, v
    return(toe, size)

size = [10,50,100,150,200,250]
f_out1 = "1electron_"
f_out2 = "points"

for sz in size:
    A, m = dots(sz)
    lam = jacobiT(A,m)
    del A
    fout = f_out1 + str(m) + f_out2
    exact = np.zeros((m))
    error = np.zeros((m))
    exact[0] = 3
    for ii in range(1,m):
        exact[ii] = exact[ii-1] + 4
        error[ii-1] = r_error(exact[ii-1],lam[ii-1])
    error[m-1] = r_error(exact[m-1],lam[m-1])
    Aux = np.hstack((lam.reshape(lam.shape+(1,)),exact.reshape(exact.shape+(1,)),error.reshape(error.shape+(1,))))
    np.savetxt(fout, Aux, fmt='%-15.6f', delimiter='   ')
    del Aux,lam, error, exact


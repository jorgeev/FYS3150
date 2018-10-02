# Buckling beam
# Libraries
import numpy as np
import math
import scipy.linalg as cpy
from jtools import *
import time
# import sys, os

def bbeam(size = 10):
    step2 = (1. / (size + 1))**2
    Aux = np.zeros(size)
    Aux[0] = 2. / step2
    Aux[1] = -1. / step2
    toe = cpy.toeplitz(Aux)
    #print(toe)
    return(toe, size)

#sz = [10,20,30,40,50,100,150,200]
sz = [10]
f_out = "results_"
title = "time_"
itere = 1 # For time the code execution

for szs in sz:
    f_name = f_out + str(szs)
    f_time = title + str(szs)
    rtime = open(f_time,"w")
    for its in range(itere):
        A,m = bbeam(szs)
        Mdiag = A[0,0]
        offD = A[0,1]
        lam,V,t = jacobiT(A,m,1)
        rtime.write(str(t)+"\n")
        # Compute exact values
        exact = np.zeros(m)
        for ii in range(m):
            exact[ii] = Mdiag + 2 * offD * math.cos((ii + 1) * math.pi * (1. / (m + 1)))
        # Create and output file with the eigenvalues
        Aux = np.hstack((lam.reshape(lam.shape+(1,)),exact.reshape(exact.shape+(1,)),abs(lam.reshape(lam.shape+(1,))-exact.reshape(exact.shape+(1,)))))
        np.savetxt(f_name, Aux, fmt='%0.6f', delimiter=',')
        del A, m, lam, V , t, exact, Aux
rtime.close()

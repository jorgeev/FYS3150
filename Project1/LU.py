#! /usr/bin/python
# -*- coding: utf-8 -*-
# LU.py

# Version: 2018.09.09.02


# Requiered libraries

import numpy as np
import scipy.linalg
from scipy.sparse import diags
import time
import sys, os
from math import exp, expm1, log10

# Solving by LU factorization using scipy library.
def LU_DSol(sz,Uy):
# Create Tridiagonal Matrix
    dd = np.array([-1 *np.ones(sz - 2),2 * np.ones(sz - 1),-1 * np.ones(sz - 2)])
    config = [-1,0,1]
    Mx = diags(dd,config).toarray()
# Applying LU solving
    t0 = time.time()
    P  = scipy.linalg.lu_factor(Mx)
    sol_LU = scipy.linalg.lu_solve(P,Uy)
    t1 = time.time()
    return(sol_LU,t1-t0)

# Equations
F = lambda x : 100.0 * exp(-10.0 * x)
exact = lambda x : 1.0 - (1 - exp(-10.0)) * x - exp(-10.0 * x)
r_err = lambda x, y : abs( (x - y) / x)

#Files output
f_out = "LU_sol104"
rtime = open("LU_time104","w")

# Iterations
it = 200

for m in range(it):
    p = 4
    sz = 10 ** p
    h = 1 / sz
    h2 = h ** 2

# Variable for functions F and exact
    x = np.linspace(0,1,sz+1)
    Ux = np.zeros((1,sz+1))
    ex = np.zeros((1,sz+1))
    ror = np.zeros((1,sz-1))

# Filling solution vector for F and exact
    for ii in range(sz+1):
        Ux[:,ii] = h2 * F(x[ii])
        ex[:,ii] = exact(x[ii])

    Ux = np.transpose(Ux[:,1:-1])
    ex = np.transpose(ex[:,1:-1])
    x = x.reshape(x.shape+(1,))
    x = x[1:-1,:]
# Run decomposition
    ss, tt = LU_DSol(sz, Ux)

# Error calculation
    for ii in range(sz-1):
        ror[:,ii] = log10(r_err(ex[ii,:],ss[ii,:]))

# Prepare output data formating
    Aux = np.hstack((x, ss, ex,np.transpose(ror)))
    np.savetxt(f_out, Aux, fmt='%0.6f', delimiter='   ')
    rtime.write(str(tt)+"\n")
    del ss, tt, ex, Aux

rtime.close()

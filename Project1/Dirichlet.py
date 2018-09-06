#! /usr/bin/python
# -*- coding: utf-8 -*-
# Dirichlet.py

# Version: 2018.09.03.03


# Requiered libraries

import numpy as np
import scipy.linalg
from scipy.sparse import diags
import time
import sys, os
from tool import *
from math import exp, expm1, log10

F = lambda x : 100.0 * exp(-10.0 * x)
exact = lambda x : 1.0 - (1 - exp(-10.0)) * x - exp(-10.0 * x)
r_err = lambda x, y : abs( (x - y) / x)

# For iterations
it = 200

f_out = "p_out100000"
rtime = open("p_time100000","w")

for m in range(it):
# Matrix & step size
    p = 5
    sz = 10 ** p
    h = 1 / sz
    h2 = h ** 2

# Variable for functions F and exact 
    x = np.linspace(0,1,sz+1)
    Ux = np.zeros((1,sz+1))
    ex = np.zeros((1,sz+1))
    ror = np.zeros((1,sz+1))

# Filling solution vector for F and exact
    for ii in range(sz+1):
        Ux[:,ii] = h2 * F(x[ii])
        ex[:,ii] = exact(x[ii])

    Uy = np.transpose(Ux[:,1:-1])

    M = 2 * np.ones((1,sz))
    S = -1 * np.ones((1,sz))
    L = S

    t0 = time.time()
    # ss = Gauss_G(M, S, L, Ux, sz)
    ss = Gauss_P(M, Ux, sz)
    # ss = LU_DSol(sz,Uy)
    t1 = time.time()
    t_f = t1 - t0
    # print(ss)
    # print(ex)
    for ii in range(1,sz):
        ror[:,ii] = log10(r_err(ex[:,ii],ss[:,ii]))
    # print(x.conj().T)
    ss = np.hstack((ss[0,:],0))
    ss = ss.reshape(ss.shape+(1,))
    x = x.reshape(x.shape+(1,))
    # print(ss)
    # print(x)
    Aux = np.hstack((x, ss,np.transpose(ex),np.transpose(ror)))
    Aux = Aux[1:-1,:]
    np.savetxt(f_out, Aux, fmt='%0.6f', delimiter='   ')
    #Aux = np.hstack((x,np.transpose(ex),np.transpose(ror)))
    # print("Iteration", m, "Runtime = ", t_f, "Seconds")
    rtime.write(str(t_f)+"\n")
    del ss, t0, t1, t_f, ex
    
rtime.close()



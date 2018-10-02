#! /usr/bin/python
# -*- coding: utf-8 -*-
# Jacobi.py

# Version: 2018.09.30.04


# Requiered libraries
# Libraries
import numpy as np
import math
import scipy.linalg as cpy
import time
from jtools import *
import sys, os

# Parameters for testing
def initialize(n = 100):
    Aux = np.zeros(n)
    Aux[0] = 2.
    Aux[1] = -1.
    toe = cpy.toeplitz(Aux)
    print(toe)
    return(toe, n)
        
A, m = initialize()
#lam,vec = jacobiT(A,m,tol=1e-10)
#print(lam)
#print("=====================================")
#print(vec)

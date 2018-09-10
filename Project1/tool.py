# Library created for FYS3150 projects
import numpy as np
import scipy.linalg
from scipy.sparse import diags

# Solution for any triangular matrix.
def Gauss_G(M, S ,L , Ux, sz):
    sol = np.zeros((1,sz))
# Forward substitution
    for ii in range(2,sz):
        M[:,ii] = M[:,ii] - (S[0,ii-1] * L[0,ii-1]) / M[:,ii-1]
        Ux[:,ii] = Ux[:,ii] - (Ux[:,ii-1]*L[0,ii-1]) / M[:,ii-1]
## 6 * (sz-3) FLOPS

# Backward substitution
    sol[:,sz-1] = Ux[:,sz-1] / M[:,sz-1]
    for ii in range(sz-2,0,-1):
        sol[:,ii] = (Ux[:,ii] - (S[:,ii] * sol[:,ii+1])) / M[:,ii]
## 3 * (sz -3) + 1 FLOPS
    return(sol)
## No. of Floating Point Operations = 9 * (sz - 3) + 1

# Solution for a tridiagonal and simetric matrix.
def Gauss_P(M, Ux, sz):
    sol = np.zeros((1,sz))
# Forward Substitution
    for ii in range(2,sz):
        M[:,ii] = M[:,ii] - 1 / M[:,ii-1]
        Ux[:,ii] = Ux[:,ii] + Ux[:,ii-1] / M[:,ii-1]
## 4 * (sz-3) FLOPS

# Backward Substitution
    sol[:,sz-1] = Ux[:,sz-1] / M[:,sz-1]
    for ii in range(sz-2,0,-1):
        sol[:,ii] = (Ux[:,ii] + sol[:,ii+1]) / M[:,ii]
    return(sol)
## 2 * (sz-3) + 1 FLOPS

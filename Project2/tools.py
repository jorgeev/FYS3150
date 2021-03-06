# Library created for FYS3150 project2. Jacobi method
import numpy as np
import math

# Search for the absolute max value, only sweep along the upper diagonal
def tracking(A,m):
        Mval = 0
        for ii in range(0,m):
            for jj in range(ii+1,m):
                if abs(A[ii,jj]) > abs(Mval):
                    Mval = A[ii,jj]
                    row = ii
                    col = jj
        return(row, col, Mval)

def rotation(A, V, m, row, col):
    # Ang calculation
    tau = (A[col,col] - A[row,row]) / (2. * A[row,col])
    t = 1. / (abs(tau) + math.sqrt(tau**2 + 1.))
    if tau < 0:
        t = -t
    c = 1. / math.sqrt(1. + t**2)
    s = c * t
    # Allocation original values
    a_row = A[row,row]
    a_col = A[col,col]
    # Rotating matrix
    A[row,row] = (a_row * c**2) + (a_col * s**2) - (2. * A[row,col] * s * c)
    A[col,col] = (a_row * s**2) + (a_col * c**2) + (2. * A[row,col] * s * c)
    A[col,row] = A[row,col] = 0.
    for ii in range(m):
        if ii != row and ii != col:
            a_row = A[ii,row]
            a_col = A[ii,col]
            A[ii,row] = A[row,ii] = a_row * c - a_col * s
            A[ii,col] = A[col,ii] = a_col * c + a_row * s

    # Calculating eigenvectors
    for ii in range(m):
        v_row = V[ii,row]
        v_col = V[ii,col]
        V[ii,row] = v_row * c - v_col * s
        V[ii,col] = v_col * c + v_row * s
    return(A,V)

def sorting(lam,V,m):
    for ii in range(m-1):
        index = ii
        val = lam[ii]
        for jj in range(ii+1,m):
            if lam[jj] < val:
                index = jj
                val = lam[jj]
        if index != ii:
            sortvalue(lam,ii,index)
            sortvector(V,ii,index)
    return(lam,V)

def sortvalue(lam,ii,jj):
    auxi = lam[ii]
    auxj = lam[jj]
    lam[ii] = auxj
    lam[jj] = auxi
    return(lam,ii,jj)

def sortvector(V,ii,jj):
    V[:,[ii,jj]] = V[:,[jj,ii]]

#m =10
#V = np.identity(m,float)
#tol = 1e-10
#mval = 1 # Just to initialize the while loop:
# Iterating to look for the a solution
#while abs(mval) > tol:
#    row,col,mval = tracking(A,m)
#    A,V = rotation(A,V,m,col,row)
     
#lam = np.diag(A)
#lam.setflags(write=1)
# Sorting results
#lam,V = sorting(lam,V,m)

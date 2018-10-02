# Library created for FYS3150 project2. Jacobi method
import numpy as np
import math
import time

def jacobiT(A,m,mode=0,tol=1e-10):
    # Search max. value
    def tracking(A,m):
        Mval = 0
        for ii in range(0,m):
            for jj in range(ii+1,m):
                if abs(A[ii,jj]) > abs(Mval):
                    Mval = A[ii,jj]
                    row = ii
                    col = jj
        return(row, col, Mval)
    
    # Calculate eigenvales and eigenvectors
    def rotation(A, m, row, col):
    
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
        A[col,row] = 0.
        A[row,col] = 0.
        for ii in range(m):
            if ii != row and ii != col:
                    a_row = A[ii,row]
                    a_col = A[ii,col]
                    A[ii,row] = a_row * c - a_col * s
                    A[row,ii] = a_row * c - a_col * s
                    A[ii,col] = a_col * c + a_row * s
                    A[col,ii] = a_col * c + a_row * s
                
        return(A)
    
    # Sort eigenvalues and eigenvectors
    def sorting(lam,m):
        for ii in range(m-1):
            index = ii
            val = lam[ii]
            for jj in range(ii+1,m):
                if lam[jj] < val:
                    index = jj
                    val = lam[jj]
            if index != ii:
                sortvalue(lam,ii,index)
        return(lam)

    def sortvalue(lam,ii,jj):
        auxi = lam[ii]
        auxj = lam[jj]
        lam[ii] = auxj
        lam[jj] = auxi
        return(lam,ii,jj)

    
    
    mval = 1
    it = 1
    if mode != 0:
        t0 = time.time()
    # Iterating to look for the a solution
    while abs(mval) > tol:
        row,col,mval = tracking(A,m)
        A = rotation(A,m,col,row)
        it += 1
        
    lam = np.diag(A)
    lam.setflags(write=1)
    if mode != 0:
        t1 = time.time()
        print("Itereations", it)
        #print("Time to compute",t1-t0)
    # Sorting results
    lam = sorting(lam,m)
    if mode == 0:
        return(lam)
    else:
        return(lam,t1-t0)

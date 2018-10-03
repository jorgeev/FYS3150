import numpy as np
import scipy.linalg as cpy
import math
from tools import tracking, rotation, sorting

def testdata1():
    A = np.array([[9.,3.,-5.,2.,-8.], [3.,9.,6.,-7.,-1.], [-5.,6.,9.,1.,9.],[2.,-7.,1.,9.,-12.],[-8.,-1.,9.,-12.,9.]])
    m = len(A)
    return(A,m)

def test_maxvalue():
    A, m = testdata1()
    row,col,maxv = tracking(A,m)
    if maxv == -12. and row == 3 and col == 4:
        return "Pass"
    else:
        return "Fail"    
    
def test_eigenvalues():
    A, m = testdata1()
    tol = 1e-10
    B = A
    V = np.identity(m,float)
    mval = 1 # Just to initialize the while loop:
    # Iterating to look for the a solution
    while abs(mval) > tol:
        row,col,mval = tracking(A,m)
        A,V = rotation(A,V,m,col,row)     
    lam = np.diag(A)
    lam.setflags(write=1)
    # Sorting results
    lam,V = sorting(lam,V,m)
    print(lam)
    # Computing values though scipy
    lam1, x1 =  cpy.eigh(B)
    tolerror = np.zeros(m)
    ii = 0
    print(lam1)
    # Comparing difference between scipy and our algorithm    
    for values in range(m):
        tolerror[values] = lam1[values] - lam[values]
        print(tolerror[values])
        if abs(tolerror[values]) <= tol:
            ii += 1

    if ii == m:
        return("Pass",tol)
    else:
        return("Fail",tol)
    
def test_ortogonality():
    A, m = testdata1()
    tol = 1e-10
    V = np.identity(m,float)
    mval = 1 # Just to initialize the while loop:
    # Iterating to look for the a solution
    while abs(mval) > tol:
        row,col,mval = tracking(A,m)
        A,V = rotation(A,V,m,col,row)     
    lam = np.diag(A)
    lam.setflags(write=1)
    # Sorting results
    lam,V = sorting(lam,V,m)
    print(V)
    v1 = V[:,0]
    v2 = V[:,1]
    v3 = V[:,2]
    v4 = V[:,3]
    v5 = V[:,4]
    ortho = np.zeros(9)
    
    o1 = np.inner(v1,v2)
    o2 = np.inner(v1,v3)
    o3 = np.inner(v1,v4)
    o4 = np.inner(v1,v5)
    o5 = np.inner(v2,v3)
    o6 = np.inner(v2,v4)
    o7 = np.inner(v2,v5)
    o8 = np.inner(v3,v4)
    o8 = np.inner(v3,v5)
    o9 = np.inner(v4,v5)
    print(o1,o2,o3,o4,o5,o6,o7,o8,o9)
    if o1 <= tol and o2 <= tol and o3 <= tol and o4 <= tol and o5 <= tol and o6 <= tol and o7 <= tol and o8 <= tol and o9 <= tol:
        return("Pass",tol)
    else:
        return("Fail",tol)

test1 = test_maxvalue()
print("Getting the max value test: ",test1)
test2,t = test_eigenvalues()
print("Difference between scipy and Jacobi algorithm with tolerance ",t,": ",test2)
test3,t = test_ortogonality()
print("Orthogonality between eigenvectors with tolerance ",t,": ",test3)


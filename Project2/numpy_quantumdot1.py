# Quantum dots
# Libraries
import numpy as np
import math
import scipy.linalg as cpy

r_error = lambda x, y : abs( (x - y) / x)

def dots(size = 5, RMin = 0.,RMax = 10.):
    # Creates a tridiagonal matrix
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
        v[ii] = r[ii]**2  # Calculating potential for the point
        toe[ii,ii] = toe[ii,ii] + v[ii] # Adding potential to the main diagonal constant
    del Aux, r, v
    return(toe, size)

#size = [500,1000,1500,2000]
size = [3000,4000,5000,10000]
f_out1 = "1electron_"
f_out2 = "points_scipy"

for sz in size:
    A, m = dots(sz)
    lam, x = cpy.eigh(A)
    
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


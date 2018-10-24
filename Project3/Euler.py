
# coding: utf-8

# In[1]:


# Libraries
import numpy as np
import math
import scipy.linalg as cpy
import time
# import time


# In[2]:


def initialize():
    size = 10000
    ftime = 10.
    step = ftime / size
    time = 0.
    # Initial values
    pi4 = 4. * math.pi**2
    x = 1.; y = 0.; z = 0.
    vx = 0.; vy = 2. * math.pi; vz = 0.
    r = math.sqrt(x**2 + y**2 + z**2)
    return(x,y,z,vx,vy,vz,r,ftime,time,step,size,pi4)
    
def to_file3d(time, x, y, z, vx, vy, vz, fout):
    aux = str('{0:0.4f},{1: 0.8f},{2: 0.8f},{3: 0.8f},{4: 0.8f},{3: 0.8f},{4: 0.8f}\n'.format(time, x, y, z, vx, vy, vz))
    fout.write(aux)


# # Euler Method

# In[3]:


# Functions for Euler Method
def Euler_method3d(x,y,z,vx,vy,vz,r,ftime,time,step,pi4,fout):
    cc="Euler_consev" #Output file for energy preservation comprovation
    ftr=open(cc,"w")
    
    
    while time <= ftime:
        x += step * vx
        y += step * vy
        z += step * vz
        vx -= step * pi4 * x / r**3
        vy -= step * pi4 * y / r**3
        vz -= step * pi4 * z / r**3
        r = math.sqrt(x**2 + y**2 + z**2)
        time += step
        to_file3d(time, x, y, z, vx, vy, vz, fout)
        # Calaculate energy and momentum preservation
        ke = kenergy(vx,vy,vz)
        pe = penergy(r)
        EE = ke + pe
        mm = momentum(x,y,z,vx,vy,vz)
        out = str('{0:0.4f},{1:0.8e},{2: 0.8e},{3: 0.8e},{4: 0.8e}\n'.format(time,ke, pe,EE, mm))
        ftr.write(out)
    fout.close()
    ftr.close()


# In[4]:


def kenergy(vx,vy,vz):#kinetic
    e_mass=5.97219e24/1988500e24
    return 0.5*e_mass*(vx**2+vy**2+vz**2)

def penergy(r):#Potential
    e_mass=5.97219e24/1988500e24
    return 4*(math.pi**2)*e_mass/r

def momentum(x,y,z,vx,vy,vz):
    e_mass=5.97219e24/1988500e24
    v3d = np.array([vx,vy,vz])
    xyz = np.array([x,y,z])
    return sum(np.cross(xyz,(e_mass*v3d))**2)
    


# In[5]:



x,y,z,vx,vy,vz,r,ftime,times,step,size,pi4 = initialize()
fnm = "Euler_"
fname = fnm + str(int(ftime)) + "years" + str(size)
fout = open(fname,"w")

to_file3d(times, x, y, z, vx, vy, vz, fout)
Euler_method3d(x,y,z,vx,vy,vz,r,ftime,times,step,pi4,fout)


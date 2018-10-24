
# coding: utf-8

# In[1]:


# Libraries
import numpy as np
import math
import scipy.linalg as cpy
import time
from astroquery.jplhorizons import Horizons
from astroquery.jplhorizons import conf
conf.horizons_server = 'https://ssd.jpl.nasa.gov/horizons_batch.cgi'

def get_planet(p_id,mss):
    obj = Horizons(id=p_id, id_type= 'majorbody', epochs={'start':'2018-10-20', 'stop':'2018-10-21','step':'1d'})
    vec = obj.vectors()
    return planet(vec['targetname'][0],np.array([vec['x'][0],vec['y'][0],vec['z'][0]]),365.25*np.array([vec['vx'][0],vec['vy'][0],vec['vz'][0]]),mss/1988500e24)


# In[2]:


class planet:
    name = "sun"
    mass = 1.
    xyz = np.array([0.,0.,0.])
    v3d = np.array([0.,0.,0.])
    acc0 = 0.
    out = open(name,"w")
    
    def __init__(self, name, xyz, velocity, mass):
        self.name = name
        self.xyz = np.array(xyz)
        self.v3d = np.array(velocity)
        self.mass = mass
        self.out = open(name,"w")
        
    def distance(planet, other_planet):
        return math.sqrt(sum((planet.xyz - other_planet.xyz)**2))

    
    def output(planet):
        aux = str('{0:0.4f},{1: 0.8f},{2: 0.8f},{3: 0.8f}\n'.format(tt, planet.xyz[0], planet.xyz[1], planet.xyz[2]))
        planet.out.write(aux)


# In[3]:


def initialsol():
    n = 1e8 # Mesh points
    ftime = 100. #Years to simulate
    h = ftime/n # Step size
    tt = 0. # initial counter
    return n, ftime, h, tt


# In[4]:


#Verlet
def accel(planets, SSystem):
    c2 = 63284.9**2
    acc = 0.
    pi4 = 4. * math.pi**2
    for ii in range(0,len(SSystem)):
        if planets.name != SSystem[ii].name:
            acc += SSystem[ii].mass / (planet.distance(planets,SSystem[ii])**3)
            if planets.name == "mercury":
                l=cpy.norm(np.cross(planets.xyz,planets.v3d))
                rc=(3*l**2)/((planet.distance(planets,SSystem[0])**2)*c2)
    acc = acc * pi4
    return acc    


# In[5]:


n, ftime, h, tt = initialsol()

#P0 = get_planet('10', 988500e24)        # Sun
#P1 = get_planet('199', 3.302e23)        # Mercury
#P2 = get_planet('299', 48.685e23)       # Venus
#P3 = get_planet('399', 5.97219e24)      # Earth
#P4 = get_planet('499', 6.4171e23)       # Mars
#P5 = get_planet('599', 1898.13e24)      # Jupiter
#P6 = get_planet('699', 5.6834e26)       # Saturn
#P7 = get_planet('799', 86.813e24)       # Uranus
#P8 = get_planet('899', 102.413e24)      # Neptun
#P9 = get_planet('999', 1.307e22)        # Pluton

P0 = planet
P1 = planet("mercury", [0.3075,0.,0.], [0.,12.44,0.], 3.003364344983656e-06)


# In[ ]:


#SSystem = [P0,P1,P2,P3,P4,P5,P6,P7,P8,P9]
#SSystem = [P0,P3,P5]
SSystem =[P0,P1]
xtr="periel"
fout=open(xtr,"w")

while tt <= ftime:
    tt += h
    for planets in SSystem:
        planets.acc0 = accel(planets, SSystem)

    for planets in SSystem:
        old_xyz = planets.xyz
        old_v3d = planets.v3d
        old_acc = -planets.acc0 * old_xyz
        for ii in range(3):
            planets.xyz[ii] = old_xyz[ii] + h*old_v3d[ii] + 0.5*h*h*old_acc[ii]
        new_acc = -planets.acc0 * planets.xyz
        for ii in range(3):
            planets.v3d[ii] = old_v3d[ii] + 0.5*h*(old_acc[ii]+new_acc[ii])
        planet.output(planets)
        if planets == "mercury":
            aux = str('{0:0.8f},{1: 0.8f},{2: 0.8f}\n'.format(planets.xyz[0],planets.xyz[1],math.atan2(planets.xyz[1]/planets.xyz[0])))
            fout.write(aux)         
for planets in SSystem:
    planets.out.close()
fout.close()


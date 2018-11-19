import numpy as np
import matplotlib.pyplot as plt

var = 1
x, y1 = np.loadtxt("Lattice_40", delimiter='   ', usecols = (0, var), unpack=True)
y2 = np.loadtxt("Lattice_60", delimiter='   ', usecols = (var), unpack=True)
y3 = np.loadtxt("Lattice_80", delimiter='   ', usecols = (var), unpack=True)
y4 = np.loadtxt("Lattice_100", delimiter='   ', usecols = (var), unpack=True)

f1 = plt.figure(1,figsize=(18, 9))
plt.plot(x,y1, '--', label='40x40')
plt.plot(x,y2, '--', label='60x60')
plt.plot(x,y3, '--', label='80x80')
plt.plot(x,y4, '--', label='100x100')
plt.tick_params(labelsize=22)
plt.xlim(2,2.4)

plt.xlabel("Temperature (T)", fontsize=22)
plt.ylabel("Energy (E/J)", fontsize=22)

plt.legend(loc='best',prop={'size': 20})
plt.grid(True)
plt.savefig("41eEnergy"+".eps", format='eps')
#plt.show()

#----------------------------------------------------------
var = 5
y1 = np.loadtxt("Lattice_40", delimiter='   ', usecols = (var), unpack=True)
y2 = np.loadtxt("Lattice_60", delimiter='   ', usecols = (var), unpack=True)
y3 = np.loadtxt("Lattice_80", delimiter='   ', usecols = (var), unpack=True)
y4 = np.loadtxt("Lattice_100", delimiter='   ', usecols = (var), unpack=True)

f2 = plt.figure(2,figsize=(18, 9));
plt.plot(x,y1, '--', label='40x40')
plt.plot(x,y2, '--', label='60x60')
plt.plot(x,y3, '--', label='80x80')
plt.plot(x,y4, '--', label='100x100')
plt.tick_params(labelsize=22)
plt.xlim(2,2.4)
plt.ylim(0,1)
plt.xlabel("Temperature (T)", fontsize=22)
plt.ylabel("Magnetisation (|M|)", fontsize=22)

plt.legend(loc='best',prop={'size': 20})
plt.grid(True)
plt.savefig("41eMagnetization"+".eps", format='eps')
#plt.show()

#----------------------------------------------------------
var = 3
y1 = np.loadtxt("Lattice_40", delimiter='   ', usecols = (var), unpack=True)
y2 = np.loadtxt("Lattice_60", delimiter='   ', usecols = (var), unpack=True)
y3 = np.loadtxt("Lattice_80", delimiter='   ', usecols = (var), unpack=True)
y4 = np.loadtxt("Lattice_100", delimiter='   ', usecols = (var), unpack=True)

f3 = plt.figure(3,figsize=(18, 9));
plt.plot(x,y1, '--', label='40x40')
plt.plot(x,y2, '--', label='60x60')
plt.plot(x,y3, '--', label='80x80')
plt.plot(x,y4, '--', label='100x100')
plt.tick_params(labelsize=22)
plt.xlim(2,2.4)
plt.xlabel("Temperature (T)", fontsize=22)
plt.ylabel("Heat Capacity (Cv/JkB)", fontsize=22)

plt.legend(loc='best',prop={'size': 20})
plt.grid(True)
plt.savefig("51eHeat_Capacity"+".eps", format='eps')
#plt.show()

#----------------------------------------------------------
var = 7
y1 = np.loadtxt("Lattice_40", delimiter='   ', usecols = (var), unpack=True)
y2 = np.loadtxt("Lattice_60", delimiter='   ', usecols = (var), unpack=True)
y3 = np.loadtxt("Lattice_80", delimiter='   ', usecols = (var), unpack=True)
y4 = np.loadtxt("Lattice_100", delimiter='   ', usecols = (var), unpack=True)

f4 = plt.figure(4,figsize=(18, 9));
plt.plot(x,y1, '--', label='40x40')
plt.plot(x,y2, '--', label='60x60')
plt.plot(x,y3, '--', label='80x80')
plt.plot(x,y4, '--', label='100x100')
plt.tick_params(labelsize=22)
plt.xlim(2,2.4)
plt.xlabel("Temperature (T)", fontsize=22)
plt.ylabel("Susceptibility (X)", fontsize=22)

plt.legend(loc='best',prop={'size': 20})
plt.grid(True)
plt.savefig("41eSusceptibility"+".eps", format='eps')
#plt.show()


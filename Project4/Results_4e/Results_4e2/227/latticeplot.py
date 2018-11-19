import numpy as np
import matplotlib.pyplot as plt

var = 1
ax0 = 2.27
axf = 2.32

x, y1 = np.loadtxt("Lattice_40", delimiter='   ', usecols = (0, var), unpack=True)
y2 = np.loadtxt("Lattice_60", delimiter='   ', usecols = (var), unpack=True)
y3 = np.loadtxt("Lattice_80", delimiter='   ', usecols = (var), unpack=True)
y4 = np.loadtxt("Lattice_100", delimiter='   ', usecols = (var), unpack=True)

f1 = plt.figure(1,figsize=(18, 9))
 
plt.plot(x,y1, '--', label='40x40')
plt.plot(x,y2, '--', label='60x60')
plt.plot(x,y3, '--', label='80x80')
plt.plot(x,y4, '--', label='100x100')
plt.xlim(ax0,axf)

plt.xlabel("Temperature (T)")
plt.ylabel("Energy (E/J)")
plt.title("Energy per Spin")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Energy"+".eps", format='eps')
plt.show()

#----------------------------------------------------------
var = 5
y1 = np.loadtxt("Lattice_40", delimiter='   ', usecols = (var), unpack=True)

y2 = np.loadtxt("Lattice_60", delimiter='   ', usecols = (var), unpack=True)
y3 = np.loadtxt("Lattice_80", delimiter='   ', usecols = (var), unpack=True)
y4 = np.loadtxt("Lattice_100", delimiter='   ', usecols = (var), unpack=True)

f2 = plt.figure(2,figsize=(5, 5));
 
plt.plot(x,y1, '--', label='40x40')
plt.plot(x,y2, '--', label='60x60')
plt.plot(x,y3, '--', label='80x80')
plt.plot(x,y4, '--', label='100x100')
plt.xlim(ax0,axf)

plt.xlabel("Temperature (T)")
plt.ylabel("Magnetization (|M|)")
plt.title("Magnetization per Spin")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Magnetization"+".eps", format='eps')
plt.show()

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
plt.xlim(ax0,axf)

plt.xlabel("Temperature (T)")
plt.ylabel("Heat Capacity (Cv/JkB)")
plt.title("Heat Capacity per Spin")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Heat_Capacity"+".eps", format='eps')
plt.show()

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
plt.xlim(ax0,axf)

plt.xlabel("Temperature (T)")
plt.ylabel("Susceptibility (X)")
plt.title("Susceptibility per Spin")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Susceptibility"+".eps", format='eps')
plt.show()


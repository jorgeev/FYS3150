import numpy as np
import matplotlib.pyplot as plt

ye1, ye2 = np.loadtxt("all_up2", delimiter='   ', usecols = (0, 1), unpack=True)
n = int(len(ye1)*0.01)
#n2 = len(ye1)
#va21 = np.var(ye1,ddof=0)
#va22 = np.var(ye2,ddof=0)

#print("Variance (Full data) =",va21, " for T = 1.0")
#print("Heat Capacity =", va21/1, " at T = 1.0")
#print("Variance (Full data) =",va22, " for T = 2.4")
#print("Heat Capacity =", va21/(2.4**2), " at T = 2.4")
ye1 = ye1[n:-1]
ye2 = ye2[n:-1]
cat = None

va21 = np.var(ye1,ddof=0)
va22 = np.var(ye2,ddof=0)
print("Variance (dropping 1%) =", va21, " for T = 1.0")
print("Variance (dropping 1%) =", va22, " for T = 2.4")
print("Average (dropping 1%) =", np.average(ye1), " for T = 1.0")
print("Average (dropping 1%) =", np.average(ye2), " for T = 2.4")

# T=1, Energy
plt.hist(ye1, cat, label='Pointing up')

plt.ylabel("Events")
plt.xlabel("Energy (E/J)")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Energy_prov10"+".eps", format='eps')
plt.show()


# T=2.4, Energy
plt.hist(ye2, cat, label='Pointing up')

plt.ylabel("Events")
plt.xlabel("Energy (E/J)")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Energy_prov24"+".eps", format='eps')
plt.show()

ye1, ye2 = np.loadtxt("all_up2", delimiter='   ', usecols = (2, 3), unpack=True)
ye1 = ye1[n:-1]
ye2 = ye2[n:-1]
cat = None

va21 = np.var(ye1,ddof=0)
va22 = np.var(ye2,ddof=0)
print("Variance (dropping 10%) =", va21, " for T = 1.0")
print("Variance (dropping 10%) =", va22, " for T = 2.4")


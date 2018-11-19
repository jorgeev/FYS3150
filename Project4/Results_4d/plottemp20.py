import numpy as np
import matplotlib.pyplot as plt

ye1, ye2, ym1, ym2 = np.loadtxt("all_up", delimiter='   ', usecols = (0, 1, 2, 3), unpack=True)
yer1, yer2, ymr1, ymr2 = np.loadtxt("random", delimiter='   ', usecols = (0, 1, 2, 3), unpack=True)

# T=1, Energy
plt.plot(ye1, '-', label='Pointing up')
plt.plot(yer1, '-', label='Random')

plt.xlabel("# Cycles")
plt.ylabel("Energy (E/J)")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Energy_T10"+".eps", format='eps')
plt.show()

# T=1, Magnetization
plt.plot(ym1, '-', label='Pointing up')
plt.plot(ymr1, '-', label='Random')

plt.xlabel("# Cycles")
plt.ylabel("Magnetization")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Magnetization_T10"+".eps", format='eps')
plt.show()

# T=2.4, Energy
plt.plot(ye2, '-', label='Pointing up')
plt.plot(yer2, '-', label='Random')

plt.xlabel("# Cycles")
plt.ylabel("Energy (E/J)")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Energy_T24"+".eps", format='eps')
plt.show()

# T=24, Magnetization
plt.plot(ym2, '-', label='Pointing up')
plt.plot(ymr2, '-', label='Random')

plt.xlabel("# Cycles")
plt.ylabel("Magnetization")

plt.legend(loc='best')
plt.grid(True)
plt.savefig("Magnetization_T24"+".eps", format='eps')
plt.show()


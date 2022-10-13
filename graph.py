import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("plotdata.dat")

plt.plot(data[:,0], data[:,1], label = "numerical")
plt.plot(data[:,0], data[:,2], label = "analitical")
plt.title("power-strain")
plt.xlabel("strain")
plt.ylabel("power")
plt.legend()
plt.grid()
plt.show()
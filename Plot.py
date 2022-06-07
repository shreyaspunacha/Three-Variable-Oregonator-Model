import numpy as np
import matplotlib.pyplot as plt

# data = np.loadtxt("Phi_0.txt")
# data = np.loadtxt("Phi_Equilibrium.txt")
# Phi_0 = data.reshape(100, 100) 
# plt.pcolormesh(Phi_0.T, vmin=0.0, vmax=1.0, cmap="jet")
# plt.colorbar()
# plt.savefig("Phi_0.png")
# plt.savefig("Phi_Equilibrium.png")

for t in range(0, 20000+500, 500):
    fname = np.loadtxt("u_%.7d.txt"%t)
    data = fname.reshape(300, 300)
    plt.pcolormesh(data.T, vmin=0, vmax=1, cmap="jet")
    plt.colorbar()
    plt.title("%0.7d"%t)
    plt.savefig("u_%0.7d.png"%t)
    plt.clf()
    fname = None
    data = None

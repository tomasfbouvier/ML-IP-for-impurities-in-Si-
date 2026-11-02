import numpy as np
import matplotlib.pyplot as plt


dmol = np.loadtxt("dmol_fit.dat")
gap = np.loadtxt("gap_fit.dat")
real = np.loadtxt("actual_gap") 



plt.plot(dmol[70:,0], dmol[70:,2], label="dmol")

plt.plot(gap[70:,0], gap[70:,2], label="gap")
plt.plot(real[:,0], real[:,1], label="ground")
plt.legend()
plt.show()

import numpy as np
import ase
from scipy.interpolate import interp1d

force = np.loadtxt("dimer_force_gap")
energy = np.loadtxt("dimer_energy_gap")

f_force = interp1d(force[:,0], force[:,1])
f_energy = interp1d(energy[:,0], energy[:,1])



x = np.linspace(0.9, 1.5, 21)
print(x, f_energy(x), f_force(x))

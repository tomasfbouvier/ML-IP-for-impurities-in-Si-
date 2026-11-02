import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from ase.calculators.abinit import Abinit
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read, write
import ase
from ase.data import atomic_numbers
import os


plt.style.use("/home/fernbouv/Desktop/mystyle.lorenzo.mplstyle")
def plot_ground_truth(file, color, idx0, derivative=0):
    
    atomss = read(file, index=f"{idx0}")
    energy = []; force = []; ds = []

    atoms = atomss[0]
    if "Si" in file:
        n=0
    else:
        n = atoms.arrays['numbers']!=14
    x0 = atoms.arrays['positions'][n]
    for i, atoms in enumerate(atomss):
        e = atoms.get_potential_energy();energy.append(e)
        #n = atoms.arrays['numbers']!=14
        f = np.linalg.norm(atoms.arrays['force'][n]); force.append(f)
        d = np.linalg.norm( atoms.arrays['positions'][n] - x0); ds.append(d)
        if "Si" in file:
            print(d)



    energy = np.array(energy)
    force  = np.array(force)
    ds = np.array(ds)

    x = np.linspace(0,9, len(energy))
    if derivative == 0:

        plt.plot(ds,energy- energy[0], 'o', mfc = "none", markeredgewidth=2, markersize=6,  color = color)
    elif derivative == 1:

        plt.plot(ds,force, 'o', mfc = "none", markeredgewidth=2, markersize=6,  color = color)
    #write("jjj", atomss, format = "extxyz")

def plot_prediction(X, color, derivative=0):

    ds, e, f = np.loadtxt(f"energy_{X}.dat", unpack=True)

    f-= f[0]
    e-=e[0]
    if derivative == 0:
        plt.plot(ds,e, color = color, linewidth = 3, label = f"{X}")
    elif derivative == 1:
        plt.plot(ds,f, color = color, linewidth = 3, label = f"{X}")


plt.figure(figsize=(6,4))

plot_prediction("P", "C0", derivative = 1)
plot_prediction("Sb", "C1", derivative = 1)
plot_prediction("Bi", "C2", derivative = 1)
#plot_prediction("Si", "C3", derivative = 0)

plot_ground_truth("P_hard.xyz", "C0", "40:", derivative = 1)
plot_ground_truth("Sb_hard.xyz", "C1", "40:", derivative = 1)
plot_ground_truth("Bi_hard.xyz", "C2", "40:", derivative = 1)
#plot_ground_truth("Si_hard.xyz", "C3", ":", derivative = 0)
plt.xlabel("Distance along given direction (Å)", fontsize=18)
#plt.ylabel("$\Delta E$ (eV)", fontsize=18)
plt.ylabel("F (eV/Å)", fontsize=18)
plt.grid(False)


handles, labels = plt.gca().get_legend_handles_labels()

# Extra handles for the data / prediction distinction
dft_handle  = Line2D([], [], color="k", marker='o', mfc = "none", markeredgewidth=2,
                     markersize=10, linestyle="None", label="DFT")
pred_handle = Line2D([], [], color="k", linestyle="-", linewidth=3,
                     label="Pred.")

#plt.legend(handles=handles + [dft_handle, pred_handle])




plt.tight_layout()
plt.show()

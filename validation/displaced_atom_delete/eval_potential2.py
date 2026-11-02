from ase.calculators.abinit import Abinit
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read, write
import ase
import numpy as np
import matplotlib.pyplot as plt
from ase.data import atomic_numbers
import os

plt.style.use("/home/fernbouv/Desktop/mystyle.lorenzo.mplstyle")


"""
parameters = {'pair_style': 'nep',
              'pair_coeff': ['* * nep.txt Si']}

files = ['nep.txt']
"""



X = "Si"


def plot_ground_truth(file, derivative=0):
    atomss = read(file, index=":40")
    energy = []
    for i, atoms in enumerate(atomss):
        e = atoms.get_potential_energy();energy.append(e)
        #f = atoms.get_forces()
        atoms.energy = e
        #atoms.force = f
    
    energy = np.array(energy)
    x = np.linspace(0,8, len(energy))
    plt.plot(x,energy- energy[0], '*', label=f'{X}')
    #write("jjj", atomss, format = "extxyz")

def calculate_one_species(X):
    print(f"Species: {X}")
    if X!="Si":
        parameters = {'pair_style': 'hybrid/overlay eam/fs tabgap nep',
                      'pair_coeff': [f'* * eam/fs ../../{X}/{X}-Si.eam.fs {X} Si',
                                     f'* * tabgap ../../{X}/{X}-Si.tabgap {X} Si no yes',
                                     '* * nep nep.txt NULL Si']} 

        files = [f'../../{X}/{X}-Si.eam.fs', f'../../{X}/{X}-Si.tabgap', 'nep.txt']

    else:
        parameters = {'pair_style': 'nep',
                      'pair_coeff': ['* * nep.txt Si']} 

        files = [f'nep.txt']

    nep = LAMMPS(files=files, **parameters)

    with open("path2.xyz") as f:
        with open("temp.xyz", "w") as f2:
            xyz = f.read()
            xyz = xyz.replace("X ", f"{X} ")
            f2.write(xyz)    

    atomss = read("temp.xyz", index=":40")
    energy = []
    for i, atoms in enumerate(atomss):
        atoms.calc = nep
        e = atoms.get_potential_energy();energy.append(e)
        f = atoms.get_forces()
        atoms.energy = e
        atoms.force = f
    
    energy = np.array(energy)[:]
    x = np.linspace(0,8, len(energy))
    plt.plot(x,energy- energy[0], '-', label=f'{X}')
    #write("jjj", atomss, format = "extxyz")

plt.figure()
for X in ["P"]:
    calculate_one_species(X)
    plot_ground_truth(f"P_hard.xyz")    
plt.xlabel("Distance travelled by the moving atom among [0.8, -0.2, 1.8]", fontsize=18)
plt.ylabel("Energy eV", fontsize=18)


plt.legend()

plt.show()  
"""

plt.figure()
force = []
for i, atoms in enumerate(atomss):
    atoms.calc = nep
    f = atoms.get_forces();force.append(np.linalg.norm(f[0]))
force = np.array(force)
x = np.linspace(0,9, len(force))
plt.plot(x,force- force[0], 'k.', label='NEP')
"""
"""
gap_force = np.loadtxt('GAP_force.csv', delimiter=',')
dft_soft_force = np.loadtxt('dft_soft_force.csv', delimiter=',')
dft_hard_force = np.loadtxt('dft_hard_force.csv', delimiter=',')
plt.plot(dft_hard_force[:,0], dft_hard_force[:,1], '.',color='r', label='DFT hard force')
plt.plot(gap_force[:,0], gap_force[:,1], '*',color='c', label='GAP force')
plt.plot(dft_soft_force[:,0], dft_soft_force[:,1], '.',color='b', label='DFT soft force')
plt.xlabel("Distance travelled by the moving atom among [0.8, -0.2, 1.8]", fontsize=18)
plt.ylabel("Force eV/Å", fontsize=18)
plt.legend()
"""
#plt.show()

os.remove("temp.xyz")

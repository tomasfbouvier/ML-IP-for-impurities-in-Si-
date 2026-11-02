from ase.calculators.abinit import Abinit
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read, write
import ase
import numpy as np
import matplotlib.pyplot as plt

parameters = {'pair_style': 'nep nep.txt',
              'pair_coeff': ['* *']}

files = ['nep.txt']

"""
pair_style hybrid/overlay eam/fs tabgap nep
pair_coeff * * eam/fs ../../${X}/${X}-Si.eam.fs Si ${X} 
pair_coeff * * tabgap ../../${X}/${X}-Si.tabgap Si ${X} no yes 
pair_coeff * * nep nep.txt Si NULL 
"""

nep = LAMMPS(parameters=parameters, files=files)




atoms = read('diamond.data', format = "lammps-data", style="atomic", Z_of_type = [0,14])
atoms.set_cell([10.9218305634, 10.9218305634,10.9218305634])
atomss = []
dxx = np.linspace(0, 5.3,  100)

direction = np.array([-1.,-1.,-1.])
direction /= np.linalg.norm(direction)

fig, axs = plt.subplots(1, 2)


with open("111_data", "a") as f:
    e0=0; f0=0
    for dx in dxx:
        positions = atoms.get_positions()
        positions[0] = dx*direction
        atoms.set_positions(positions)
        atoms.set_calculator(nep)
        e = atoms.get_potential_energy()

        force = np.linalg.norm(atoms.get_forces()[0])
        if dx==0:
            e0=0 #e
            f0=force

        axs[0].plot(dx, e-e0, 'b.')
        axs[1].plot(dx, force-f0, 'b.')
    """
    for dx in dxx:
        positions = atoms.get_positions()
        positions[0] = dx*direction
        atoms.set_positions(positions)

        atoms.set_calculator(nep_old)
        e = atoms.get_potential_energy()
        force = np.linalg.norm(atoms.get_forces()[0])
        if dx==0:
            e0=0#e
            f0=force



        axs[0].plot(dx, e-e0, 'r.')

        axs[1].plot(dx, force-f0, 'r.')
    """ 
        write(f, atoms, format="xyz")
    f.close()
plt.show()

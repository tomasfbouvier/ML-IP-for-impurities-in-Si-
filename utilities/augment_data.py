import os
from ase.io.vasp import read_vasp
from ase.io import read, write
from ase.calculators.vasp import Vasp
import numpy as np
import ase.neighborlist as nl
from ase.neighborlist import get_connectivity_matrix
from ase import Atoms
from ovito.io.ase import ase_to_ovito
from ovito.data import NearestNeighborFinder
import ovito 
from ovito.modifiers import *
"""
for j in range(N):
    atoms = atoms0.copy()

    atoms[0].position[0] = atoms[0].position[0] + 5.43/(1.3*N)*j
    atoms[0].symbol = 'C'
    
    atoms.set_cell(atoms0.get_cell()*np.random.uniform(0.99,1.01) , scale_atoms=True)
    atoms.rattle(stdev=0.05, seed = 723762+j)


    write_extxyz('C_tde.xyz', atoms, append=True)

atoms0 = ase.io.read('Si.xyz' )
N=10
for j in range(-int(3*N/2),int(N/2)):
    atoms = atoms0.copy()

    atoms[0].position[0] = atoms[0].position[0] + 5.43/(1.4*np.sqrt(3)*N)*j
    atoms[0].position[1] = atoms[0].position[1] + 5.43/(1.4*np.sqrt(3)*N)*j
    atoms[0].position[2] = atoms[0].position[2] + 5.43/(1.4*np.sqrt(3)*N)*j

    atoms[0].symbol = 'C'
    
    atoms.set_cell(atoms0.get_cell()*np.random.uniform(0.99,1.01) , scale_atoms=True)
    atoms.rattle(stdev=0.05, seed = 723762+j)


    write_extxyz('C_tde.xyz', atoms, append=True)
"""
#print(atoms0)
#print(nl.natural_cutoffs(atoms0))
for i in range(10):
    #atoms0 = read('dumbell.xyz', index=0)
    #cell = atoms0.get_cell()
    atoms=read('weird3.xyz', index=0)
    scaling_factor = 1.0 #0.9 + (1.1 - 0.9) * (i / 49)  # 49 to ensure 50 steps are evenly spaced
    #atoms.set_cell(atoms0.get_cell() * scaling_factor, scale_atoms=True)
    atoms.rattle(stdev=0.14, seed=2723451*(1+i))
    write('augmented2.xyz', atoms, append=True, format="extxyz")



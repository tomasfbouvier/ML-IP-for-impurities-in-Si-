import os
from ase.io.vasp import read_vasp
from ase.io import read, write
from ase.calculators.vasp import Vasp
import numpy as np
import ase.neighborlist as nl
from ase.neighborlist import get_connectivity_matrix
from ase import Atoms, Atom
from ovito.io.ase import ase_to_ovito
from ovito.data import NearestNeighborFinder, CutoffNeighborFinder
import ovito 
from ovito.modifiers import *
for i in range(50):
    j=0    
    while True:
        atoms0=read('Sb_subs.xyz', index=0)
        atoms = atoms0.copy()
        v = np.random.rand(3)
        v/= np.linalg.norm(v)
        v*= 5*np.random.rand()
        atoms.append(Atom('Si', v ))
        #atoms[0].position[0] = atoms[0].position[0] + 5.43/(1.3*N)*j
        #atoms[0].symbol = 'Si'
        
        #atoms.set_cell(atoms0.get_cell()*np.random.uniform(0.99,1.01) , scale_atoms=True)


        atoms.rattle(stdev=0.07, seed=723762*(1+i)*(1+j))
        data = ase_to_ovito(atoms)
        finder = CutoffNeighborFinder(2.2, data).find_all()[1]

        sb_idx = np.where(data.particles['Particle Type'][...]==2)[0][0]
        finder2 = CutoffNeighborFinder(2, data)
        gen = finder2.find(sb_idx)
        count = sum(1 for _ in gen)
        #if count == 0:
        #    write('augmented.xyz', atoms, append=True, format="extxyz")
        #    break


        j+=1
        if np.min(np.linalg.norm(finder, axis=1)) > 1.7 and count == 0:
            break
        else:
            print(i)
        del finder
    write('augmented.xyz', atoms, append=True, format="extxyz")



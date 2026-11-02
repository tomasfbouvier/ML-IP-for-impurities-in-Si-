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

atoms0=read('../train.xyz', index="706:756")

for atoms in atoms0: 

    data = ase_to_ovito(atoms)
    sb_idx = np.where(data.particles['Particle Type'][...]==2)[0][0]
    finder = CutoffNeighborFinder(2, data)
    gen = finder.find(sb_idx)
    count = sum(1 for _ in gen)
    if count == 0:
        write('augmented2.xyz', atoms, append=True, format="extxyz")



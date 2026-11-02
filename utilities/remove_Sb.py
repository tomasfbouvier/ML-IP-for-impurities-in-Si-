import argparse
from ase.io import read, write
from ase.calculators.lammpsrun import LAMMPS
from ase import Atoms
from ovito.io.ase import ase_to_ovito
import numpy as np
import os


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("database", type=str)
    args = parser.parse_args()

    db = read(args.database, index=':', format="extxyz")

    for atoms in db:   
        atoms2 = atoms[[atom.index for atom in atoms if atom.symbol=='Si']]        
        write(f"{args.database}_removed", atoms2, format="extxyz", append = True)

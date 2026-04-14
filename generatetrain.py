import os
import numpy as np
import string
import subprocess
import argparse
from ase.io import read, write
from ase import Atoms
import glob



def main(args):
    try:
        os.remove("train.xyz")
    except: pass

    basedir = 'databases/X_removed/db.xyz'
   
    with open(basedir, 'r') as file:
      filedata = file.read()

    filedata = filedata.replace('energy', 'Energy')

    with open(basedir, 'w') as file:
      file.write(filedata)

    file_path = os.path.join(args.directory, "db.xyz")
   
    with open(file_path, 'r') as file:
      filedata = file.read()

    filedata = filedata.replace('energy', 'Energy')

    with open(file_path, 'w') as file:
      file.write(filedata)

    db = read(file_path, index=':', format="extxyz")
    db_base = read(basedir, index=':', format="extxyz")

    #There is something very stupid happening with the naming of metadata. "energy" is not read by ASE's file reader but "Energy" is
    #TODO: Find a way to replace "energy" by "Energy" upon reading the file and undo the change upon writting 

    for i in range(len(db)):

        atoms = db[i]
        atoms1= db_base[i]
        impurities = atoms[[atom.index for atom in atoms if atom.symbol!='Si']]
        
        Ni = len(atoms[[atom.index for atom in atoms if atom.symbol!='Si']])
        N = len(atoms)-Ni

        del atoms[[atom.index for atom in atoms if atom.symbol!='Si']]
        
        pos = atoms.positions; symb = atoms.symbols; cell = atoms.cell
        config_type  = atoms.info['config_type']

        atoms2 = Atoms(positions = pos, symbols = symb, pbc = True, cell = cell)
        atoms2.info['config_type'] = config_type 
        
        try:
            e = atoms.info['Energy'] - atoms1.info['Energy'] #- N*e0
            atoms2.info['energy'] = e
        
        except:
            raise ValueError(f"{i}, I am missing energy!")

        
        try:
            f = atoms.get_array('force') - atoms1.get_array('force')
            atoms2.new_array('force', f)
        
        except:
            pass
        
        try:
            virial = atoms.info['virial']-atoms1.info['virial'] #- N*virial0
            atoms2.info['virial'] = virial
        
        except:
            pass
       
        atoms2 +=impurities
        
        write(f"train.xyz", atoms2, format="extxyz", append = True)

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', required=True)

if __name__=="__main__":
    args = parser.parse_args()
    main(args)

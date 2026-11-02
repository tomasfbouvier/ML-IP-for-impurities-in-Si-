from ase.calculators.abinit import Abinit
from ase.calculators.lammpsrun import LAMMPS
from ase.io import read, write
import ase
import numpy as np
import matplotlib.pyplot as plt
from ase.data import atomic_numbers
import os


def calculate_one_species(X="P", derivative = 1):
    
    if X!="Si":

        parameters = {'pair_style': 'hybrid/overlay eam/fs tabgap nep',
                      'pair_coeff': [f'* * eam/fs ../../{X}/{X}-Si.eam.fs {X} Si',
                                     f'* * tabgap ../../{X}/{X}-Si.tabgap {X} Si no yes',
                                     '* * nep nep.txt NULL Si']} 

        files = [f'../../{X}/{X}-Si.eam.fs', f'../../{X}/{X}-Si.tabgap', 'nep.txt']

     
    else: 
        parameters = {'pair_style': 'nep',
                      'pair_coeff': ['* * nep.txt Si']} 

        files = ['nep.txt']

    
    calc = LAMMPS(files=files, **parameters)
    with open("Sb_subs.data") as f:
        with open("temp.data", "w") as f2:
            xyz = f.read()
            xyz = xyz.replace("X", f"{X} ")
            f2.write(xyz)   

    if X == "Si":
        with open("dia.data") as f:
            with open("temp.data", "w") as f2:
                xyz = f.read()
                f2.write(xyz)   



    atoms = read("temp.data", format = "lammps-data")
    atoms.calc = calc

    dmax = 8.5
    direction = np.array([-0.0839203652, 0.0209800913, 0.1888208218])
    direction = direction/np.linalg.norm(direction)
    
    es = []
    print(direction)
    for d in np.linspace(0, dmax, 100):
        atoms.positions[0] = d*direction
        print(atoms[0], atoms.get_potential_energy())
        
        es.append(atoms.get_potential_energy())



        write(f"predict_{X}.xyz", atoms,  append = True) 

    plt.figure()
    plt.plot(es)

    plt.show()
    """


    atoms.calc = calc

    

    energy = []; force = []

    for i, atoms in enumerate(atomss):
        atoms.calc = nep
        e = atoms.get_potential_energy();energy.append(e)
        fs = atoms.get_forces()
        atoms.energy = e
        n = atoms.arrays['numbers']!=14
        f = np.linalg.norm(fs[n]); force.append(f) 
    energy = np.array(energy)
    force = np.array(force)
    
    if derivative == 0 :
        x = np.linspace(0,9, len(energy))
        y = energy - energy[0]
    elif derivative == 1:
        x = np.linspace(0,9, len(force))
        y = force
        print(y)
    from scipy.interpolate import make_interp_spline

    x_smooth = np.linspace(x.min(), x.max(), 300)
    y_smooth = make_interp_spline(x, y, k=1)(x_smooth)

    line, = plt.plot(x_smooth, y_smooth, '-', label=f'{X}', linewidth=3)
    color = line.get_color()
    
    return color

    """


calculate_one_species("Si")

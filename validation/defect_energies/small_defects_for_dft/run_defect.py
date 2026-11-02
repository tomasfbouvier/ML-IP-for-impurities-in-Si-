import os
import subprocess
from ase.data import atomic_masses, atomic_numbers


def parse_energy_from_log():

    import re

    e_pair = None

    with open("log.lammps", "r") as f:
        for line in f:
            if re.match(r'ec', line):
                parts = line.split()
                if len(parts) >= 3:
                    e_pair = float(parts[2])  # E_pair is 3rd column
                    N = float(parts[3])
    return e_pair, N



def calculate_formation_energies(impurity):
    mydic = {}
    c=0
    mass = atomic_masses[atomic_numbers[impurity]]
    for file in os.listdir("."):
        if file.endswith(".data"):
            print(file)
            subprocess.run(f"mpirun -np 6 lmp -in in.relax -v input {file} -v X {impurity} -v mass {mass} -v output {file}.out",
                           shell=True,  stdout=subprocess.DEVNULL)
            
            mydic[os.path.splitext(file)[0]] = parse_energy_from_log() 

    formation_energies = {}
    for key in mydic.keys():
        if key in ["diamond"]: continue
        if "Sb_" not in key: continue
        e_pair, N = mydic[key] 
        e_pair_base, N_base = mydic[key.replace("Sb_", "")]
        #formation_energies[key] = e_pair -mydic["Sb_subs"][0]-1/216*mydic["diamond"][0]*(N-216)
        formation_energies[key] = e_pair - e_pair_base - mydic["Sb_subs"][0] + mydic["diamond"][0]

    return formation_energies


for impurity in ["P", "Sb", "Bi"]:
    print(impurity, calculate_formation_energies(impurity))

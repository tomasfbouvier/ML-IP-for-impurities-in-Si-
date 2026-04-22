import os
import subprocess
from ase.data import atomic_masses, atomic_numbers
from adapt_data import adapt_data
from parse_energy_barriers import parse_barrier

def relax(impurity, file_in, file_out):

    mass = atomic_masses[atomic_numbers[impurity]]
    subprocess.run(f"mpirun -np 6 lmp -in in.relax2 -v input {file_in} -v X {impurity} -v mass {mass} -v output {file_in}.out",
                   shell=True,  stdout=subprocess.DEVNULL)
    
    subprocess.run(f"mpirun -np 6 lmp -in in.relax2 -v input {file_out} -v X {impurity} -v mass {mass} -v output {file_out}.out",
                   shell=True,  stdout=subprocess.DEVNULL)
        
    adapt_data(f"{file_out}.out")    

    return


def neb(impurity, initial, final):
    
    mass = atomic_masses[atomic_numbers[impurity]]
    npart = 5
    subprocess.run(f"mpirun -np {npart} --oversubscribe lmp -partition {npart}x1 -v dump_file kkk -v num_replicas {npart} -in in.neb -v X {impurity} -v mass {mass} -v initial {initial} -v final {final} > log.neb",
                   shell=True,  stdout=subprocess.DEVNULL)
    
    parse_barrier()
    
    return

if __name__=="__main__":
    initial = "initial.data"
    final = "final2.data"
    
    relax("P", initial, final)
    
    neb("P", "initial.data.out", "final.data")

    #for impurity in ["P", "Sb", "Bi"]:
    #    print(impurity, calculate_formation_energies(impurity))














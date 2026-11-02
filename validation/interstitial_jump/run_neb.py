import os
import subprocess
from ase.data import atomic_masses, atomic_numbers
from adapt_data import adapt_data
from parse_energy_barriers import parse_barrier
import matplotlib as mpl


def relax(impurity, file_in, file_out):

    mass = atomic_masses[atomic_numbers[impurity]]
    subprocess.run(f"mpirun -np 6 lmp -in in.relax2 -v input {file_in} -v X {impurity} -v mass {mass} -v output {file_in}.out",
                   shell=True,  stdout=subprocess.DEVNULL)
    
    subprocess.run(f"mpirun -np 6 lmp -in in.relax2 -v input {file_out} -v X {impurity} -v mass {mass} -v output {file_out}.out",
                   shell=True,  stdout=subprocess.DEVNULL)
        
    subprocess.run(f"cp {file_in}.out initial.data", shell=True, stdout=subprocess.DEVNULL)
    adapt_data(f"{file_out}.out", refference= f"{file_in}.out")    
    return


def neb(impurity, initial, final, npart=5):
    
    mass = atomic_masses[atomic_numbers[impurity]]
    subprocess.run(f"mpirun -np {npart} --oversubscribe lmp -partition {npart}x1 -v dump_file kkk -v num_replicas {npart} -in in.neb -v X {impurity} -v mass {mass} -v initial {initial} -v final final.data > log.neb",
                   shell=True,  stdout=subprocess.DEVNULL)
        
    x,y = parse_barrier(plot = False)
    
    return x, y

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline


plt.style.use("/home/fernbouv/Desktop/mystyle.lorenzo.mplstyle")
"""
if __name__=="__main__":
    initial = "2NN.data_bis"
    final = "3NN.data_bis"
    relax("P", initial, final)
    x, y = neb("P", f"{initial}.out", "final.data", 5)
"""

if __name__ == "__main__":

    states = ["0NN.data", "1NN.data", "2NN.data", "3NN.data"]
    


    impurities = ["P","Sb", "Bi"]


    plt.figure(figsize = (6,4))
    markers = ["s", "D", "^"]
    markers = ["s-", "D-", "^-"]
    colors = ["C0", "C1", "C2"]

    
    points = [0.0, 1.0, 2.0, 4.0]
    plt.xticks(points,
    ["A", "B", "C", "D"])

    for i in range(len(impurities)):

        all_x = []
        all_y = []

        x_offset = 0.0
        impurity = impurities[i]
        marker = markers[i]
        color = colors[i]

        min_y = 0
        for j in range(1, len(states)):
            initial = states[j - 1]
            final = states[j]
            
            initial2 = initial.removesuffix(".data")
            final2 = final.removesuffix(".data")

            fname = f"{impurity}_{initial2}_{final2}.txt"

            # ---- Skip calculation if file exists ----
            if not os.path.exists(fname):
                relax(impurity, initial, final)
                
                if "3NN" in final2: n = 17
                else: n= 17
                x, y = neb(impurity, f"{initial}.out", "final.data", npart = n)

                with open(fname, "w") as f:
                    for xi, yi in zip(x, y):
                        f.write(f"{xi}, {yi}\n")
            else:
                # ---- Read existing data ----
                x, y = [], []
                with open(fname, "r") as f:
                    for line in f:
                        xi, yi = map(float, line.strip().split(","))
                        x.append(xi)
                        y.append(yi)

            x = np.array(x)
            y = np.array(y)
            
            x = np.linspace(points[j-1], points[j], len(x))
            y = np.array(y)
           
            if j == 1: min_y = min(y) 
            #min_y = min(0, min(y))

            y-=min_y
            plt.plot(x, y, markers[i], label=f"{impurity}", color = colors[i], linewidth=3)
            # ---- Shift x so segments connect ----
            x_shifted = x + x_offset

            # update offset for next segment
            x_offset = x_shifted[-1]

            all_x.append(x_shifted)
            all_y.append(y)

        # ---- Concatenate all segments ----
        all_x = np.concatenate(all_x)
        all_y = np.concatenate(all_y)
        
        all_x, unique_indices = np.unique(all_x, return_index=True)
        all_y = all_y[unique_indices]
        all_y[:]-=all_y[0]
        
        # ---- Cubic spline ----
        #cs = make_interp_spline(all_x, all_y, 1)
        #x_fine = np.linspace(all_x.min(), all_x.max(), 500)
        #y_fine = cs(x_fine)
        
        print(impurity)
        #all_x = [0.,   0.25, 0.5,  0.75, 1. , 1.25, 1.5, 1.75, 2.0]
        print(all_x, all_y)
        all_x = np.linspace(0, len(all_y), len(all_y)) 
        #plt.plot(all_x, all_y, marker, label=f"{impurity}", c = color)
        #plt.plot(x_fine, y_fine, "-", c = color)

    x = np.array([0.0, 0.75, 1.0, 1.25, 2.0, 3.0, 4.0])
    P = np.array([-1173.87133525, -1173.88941904, -1173.94525876, -1173.88941904, -1173.87133525, -1173.87133525,-1173.87133525])
    P -= P[0]
    
    Sb = np.array([-1171.20017128, -1169.35530386, -1169.97592997, -1169.27117777, -1171.25066641, -1170.69875198, -1171.27993963])
    Sb -= Sb[0]
    Bi = np.array([-1170.2776458, -1168.58565725, -1168.64329067, -1168.65746428, -1170.31712972,-1169.94805559, -1170.31943924])
    Bi -= Bi[0]

    ### END OUR DFT DATA###

    plt.plot(x, P, "s--", mfc='none', markersize = 7,  color="black", label="OUR")   
    plt.plot(x, Sb, "D--", mfc='none', markersize = 7,  color="black", label="OUR")   
    plt.plot(x, Bi, "^--", mfc='none', markersize = 7,  color="black", label="OUR")   


    plt.xlabel("Reaction coordinate")
    plt.ylabel("Energy (eV)")
    plt.legend()

    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # later duplicates overwrite earlier ones
    plt.legend(by_label.values(), by_label.keys(), labelspacing=0.1, columnspacing=0.5, handletextpad=0.1 )


    plt.show()
"""
if __name__=="__main__":
    

    states = ["0NN.data", "1NN.data"]
    
    for i in range(1,2):
        
        initial = states[i-1]
        final = states[i]
        impurity = "P"    
        relax(impurity, initial, final)

        x, y = neb(impurity, f"{initial}.out", "final.data")
        
        initial = initial.removesuffix(".data").removesuffix(".out")
        final = final.removesuffix(".data").removesuffix(".out")
        with open(f"{impurity}_{initial}_{final}.txt", "w") as f:
            for xi, yi in zip(x, y):
                f.write(f"{xi}, {yi}\n")

            f.close()
    
    plt.figure()
    """

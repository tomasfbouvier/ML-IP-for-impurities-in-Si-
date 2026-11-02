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
        
    adapt_data(f"{file_out}.out")    

    return


def neb(impurity, initial, final, npart=5):
    
    mass = atomic_masses[atomic_numbers[impurity]]
    print(f"mpirun -np {npart} --oversubscribe lmp -partition {npart}x1 -v dump_file kkk -v num_replicas {npart} -in in.neb -v X {impurity} -v mass {mass} -v initial {initial} -v final {final} > log.neb")
    subprocess.run(f"mpirun -np {npart} --oversubscribe lmp -partition {npart}x1 -v dump_file kkk -v num_replicas {npart} -in in.neb -v X {impurity} -v mass {mass} -v initial {initial} -v final {final} > log.neb",
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
    initial = "0NN.data"
    final = "1NN.data"
    relax("Bi", initial, final)
    x, y = neb("Bi", f"{initial}.out", "final.data", 7)
"""
if __name__ == "__main__":

    states = ["0NN.data", "1NN.data", "3NN.data", "4NN.data"]
    


    impurities = ["P", "Sb", "Bi"]

    plt.figure(figsize = (6,4))
    markers = ["s-", "D-", "^-"]
    #colors = ["red", "green", "blue"]
    colors = ["C0", "C1", "C2"]

    
    points = [0.0, 1.0, 3.0, 4.0]
    plt.xticks(points,
    ["A", "A'", "C", "D"])


    for i in range(len(impurities)):

        all_x = []
        all_y = []

        x_offset = 0.0
        impurity = impurities[i]
        marker = markers[i]
        color = colors[i]
        min_y = 0.0 
        for j in range(1, len(states)):
            initial = states[j - 1]
            final = states[j]
            print(points[j-1], points[j])


            initial2 = initial.removesuffix(".data")
            final2 = final.removesuffix(".data")

            fname = f"{impurity}_{initial2}_{final2}.txt"

            # ---- Skip calculation if file exists ----
            if not os.path.exists(fname):
                relax(impurity, initial, final)
                
                if "0NN" in initial2: n = 17
                else: n=17
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
            x = np.linspace(points[j-1], points[j], len(x))
            y = np.array(y)
            if min(y)<min_y: min_y=min(y)

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
        all_y[:]-=min(all_y)
        
        # ---- Cubic spline ----
        cs = make_interp_spline(all_x, all_y, 3)
        x_fine = np.linspace(all_x.min(), all_x.max(), 500)
        y_fine = cs(x_fine)
        
        print(impurity)
        #all_x = [0.,   0.25, 0.5,  0.75, 1.,   1.25, 1.5,  1.75, 2.  ]
        #print(all_x, all_y)
        #plt.plot(all_x, all_y, marker, c = color)
        #plt.plot(x_fine, y_fine, "-", label=f"{impurity}", c = color, linewidth=3)

    
    plt.xlabel("Reaction coordinate")
    plt.ylabel("Energy (eV)")
    
    #"""

    ##############THEORETICAL VALUES###############
    x = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0])
    P = np.array([0.0, 0.70833, 1.345898, 0.700, 0.0 , 0.8355222340301367, 0.7053983328316675, 0.8473717834524023, 0.8])
    print(x)
    #x = np.linspace(all_x.min(), all_x.max(), len(P))
    x = np.array([0., 0.25, 0.5,  0.75, 1.0,   1.5, 2.,  2.5, 3.   ])
    plt.plot(x,P, "s-",c="black", markersize = 7)
    
    ##############SPECIAL CASE: Si self vacancy diffusion#####################

    Si = np.array([-1162.6083, -1162.4903,  -1162.4036, -1162.49, -1162.6083, -1162.4036, -1162.6083, -1162.4036, -1162.6083])
    
    Si-= Si[0]

    #x = np.arange(len(Si))/4
    plt.plot(x,Si, "o-",c="C3", label="Si")
    
    Sb = np.array([0.0, 0.154, 0.249, 0.150, 0.0, 0.944, 0.903, 1.083, 0.992])

    plt.plot(x,Sb, "D-",c="black", label = "DFT")
    
    Si = np.array([0.0,0.2, 0.3, 0.2, 0.0, 0.3,0.0, 0.3, 0.0])
    plt.plot(x,Si, "o-",c="black", markersize= 7)
    #plt.plot(x,Si, "--",c="black", label = "DFT")
    ##############END OF THEORETICAL VALUES ##############################
    #"""
    

    ### OUR DFT DATA###
    

    x = np.array([0.0, 0.5, 1.0, 2.5, 3.0, 3.5, 4.0])

    Bi = np.array([-1161.19, -1161.34, -1161.18, -1159.76, -1159.61, -1159.50, -1159.53])
    Bi-=Bi[0]
    Sb = np.array([-1161.64, -1161.67, -1161.64, -1160.71, -1160.57964122, -1160.37, -1160.51])
    Sb-=Sb[0]
    P = np.array([-1163.33264132, -1162.19756723,-1163.33392446, -1162.86021026, -1162.61115714, -1162.48749253,-1162.55255395 ])
    P -= P[0]

    ### END OUR DFT DATA###

    plt.plot(x, P, "s--", mfc='none', markersize = 7,  color="black", label="OUR")
    plt.plot(x, Sb, "D--", mfc='none', markersize = 7,  color="black", label="OUR")
    plt.plot(x, Bi, "^--", mfc= 'none', markersize = 7, color = "black", label="OUR")
    plt.ylim(0, 2.0)
    plt.legend(ncols=3, labelspacing=0.1, columnspacing=0.5, handletextpad=0.1)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # later duplicates overwrite earlier ones
    plt.legend(by_label.values(), by_label.keys(),ncols=3, labelspacing=0.1, columnspacing=0.5, handletextpad=0.1 )
    plt.grid(False)
    plt.savefig("energy_barriers.pdf")
    plt.show()


"""
if __name__=="__main__":
    

    states = ["0NN", "1NN", "3NN"]
    
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


    #for impurity in ["P", "Sb", "Bi"]:
    #    print(impurity, calculate_formation_energies(impurity))














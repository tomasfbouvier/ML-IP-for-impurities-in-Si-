from collections import Counter
import numpy as np
import os
import matplotlib.pyplot as plt
from ovito.io import import_file
from ovito.modifiers import *

plt.style.use("/home/fernbouv/Desktop/mystyle.lorenzo.mplstyle")
def plot_one_species(species, offset, width=0.15):

    pipeline = import_file('/home/fernbouv/Desktop/ML_IP_impurities/validation/reconstruction/P/final.0',
                           atom_style='atomic')
    pipeline.modifiers.append(IdentifyDiamondModifier())

    states = []
    for root, dirs, files in os.walk('.'):
        for file in files:
            path = os.path.join(root, file)
            if species in path and "final" in path:
                pipeline.source.load(path)
                data = pipeline.compute()
                defect = data.particles["Structure Type"][data.particles["Particle Type"] == 2][0]

                states.append("SUBSTITUTIONAL" if defect == 1 else "DEFECT")

    counts = Counter(states)

    labels = ["SUBSTITUTIONAL", "DEFECT"]
    x = np.arange(len(labels))

    plt.bar(x + offset,
            [counts.get(l, 0) for l in labels],
            width=width,
            label=species)

plt.figure()

plot_one_species("P", -0.1)
plot_one_species("Sb", 0.1)
plot_one_species("Bi", 0.3)

plt.xticks(np.arange(2), ["SUBSTITUTIONAL", "DEFECT"])
plt.legend(title="Atomic Species")
plt.xlabel("State")
plt.ylabel("Percentage")
#plt.show()
plt.savefig("reconstruction.png")

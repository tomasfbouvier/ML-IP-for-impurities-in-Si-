import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

fname = sys.argv[1] if len(sys.argv) > 1 else "energy_map.txt"

#plt.style.use("/home/fernbouv/Desktop/mystyle.lorenzo.mplstyle")
# Read data (lines starting with '#' are skipped) and E0 from the header
data = np.loadtxt(fname)
dx, dy, e = data[:, 0], data[:, 1], data[:, 2]

e0 = None
with open(fname) as f:
    for line in f:
        if line.startswith("# E0"):
            e0 = float(line.split("=")[1].split()[0])
            break

# Energy relative to the undisplaced configuration (in meV)
if e0 is None:
    e0 = e.min()
z = (e - e0) * 1000.0

# Reshape onto the regular grid
xs, ys = np.unique(dx), np.unique(dy)
X, Y = np.meshgrid(xs, ys)
Z = np.full(X.shape, np.nan)
ix = np.searchsorted(xs, dx)
iy = np.searchsorted(ys, dy)
Z[iy, ix] = z

fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection="3d")
surf = ax.plot_surface(X, Y, Z, cmap=cm.viridis, edgecolor="k",
                       linewidth=0.2, antialiased=True)
ax.set_xlabel("dx (Å)")
ax.set_ylabel("dy (Å)")
ax.set_zlabel("ΔE (meV)")
#ax.set_title("Energy landscape for atom 109")
fig.colorbar(surf, shrink=0.6, pad=0.1, label="ΔE (meV)")
plt.rcParams["axes3d.mouserotationstyle"] = "azel"   # matplotlib >= 3.10
plt.tight_layout()
#plt.savefig("energy_map_3d.png", dpi=200)
plt.show()

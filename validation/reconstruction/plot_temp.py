import matplotlib.pyplot as plt
import numpy as np



plt.style.use("/home/fernbouv/Desktop/mystyle.lorenzo.mplstyle")
time, temp = np.loadtxt("temperature_profile.txt", unpack=True)

plt.figure(figsize=(6,3))
plt.plot(time, temp, c="blue")
plt.xlabel("Time (ns)")
plt.ylabel("Temperature (K)")
plt.tight_layout()
plt.savefig("temperature_profile.pdf")

##This code plot the 2D free energy surface with distance between two residues and fraction of native contacts as order parameters##
##The distance data and the fraction of native contacts can be calculated using a different script##
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from MDAnalysis.analysis import distances, rms, polymer
import sys

a = sys.argv[1]
b = sys.argv[2]

Q_frac = np.genfromtxt(a)
dist = np.genfromtxt(b)

k_B = 0.0019872041 
temp = 300  

hist, x, y = np.histogram2d(Q_frac, dist, bins=100, density=True)
free_energy = -k_B * temp * np.log(hist + 1e-10)  

free_energy = np.ma.masked_invalid(free_energy)
X, Y = np.meshgrid(x[:-1], y[:-1])

plt.contourf(X, Y, free_energy, levels=20, cmap='viridis')
plt.colorbar(label='Free Energy (kcal/mol)')
#plt.xlabel('Q_frac')
#plt.ylabel('dist')
#plt.ylim(0.0,60.0)
#plt.title('2D Free Energy Surface')
plt.show()

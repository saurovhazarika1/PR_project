import numpy as np
import sys
import pyemma as pym
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.axis import Axis

a = sys.argv[1]
b = sys.argv[2]

data1 = np.genfromtxt(a)
data2 = np.genfromtxt(b)

RMSD1 = np.array(data1)
RMSD2 = np.array(data2)
#rg = data2[:,1]

#fig, ax = plt.subplots()

pym.plots.plot_free_energy(RMSD1, RMSD2)
#plt.xlim(50,160)
#plt.ylim(-5,5)
#plt.xlabel('RMSD1')
#plt.ylabel('RMSD2')
#plt.ylim((0, 16))
#plt.xlim((0,40))
#ax.xaxis.set_minor_locator(tick.AutoMinorLocator())
#plt.title('log(P) plot of 4a2j')
#plt.savefig('log_prob_anta.png')
plt.show()

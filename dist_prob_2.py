import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import seaborn as sns

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob

import sys

trajectory_files = glob('4a2j/4a2j_strip_*.dcd')
print(trajectory_files)
topology_file = '4a2j/4a2j_strip.prmtop'


#sel_basic = "(resname ARG LYS) and (name NH1 NE NH2 NZ)"
#sel_acedic = "(resname ASP GLU) and (name OD1 OD2 OE1 OE2)"

sel_basic_1 = "(resid 217) and (name CA)"
sel_acedic_1 = "(resid 41) and (name CA)"

sel_basic_2 = "(resid 42) and (name NH2)"
sel_acedic_2 = "(resid 222) and (name OE2)"

 
timeseries = []

for traj in trajectory_files:
    u = mda.Universe(topology_file, traj)
    for ts in u.trajectory:
        acedic = u.select_atoms(sel_acedic_2)
        basic = u.select_atoms(sel_basic_2)
        dist = contacts.distance_array(acedic.positions, basic.positions)
        print(np.array(dist).shape)
        #dist_f = np.linalg.norm(dist[0]-dist[1])
        #print(dist_f)
        timeseries.append(dist[0][0])
        #n_contacts = contacts.Contacts(u, select=(sel_acedic_1, sel_basic_1))
        #print(n_contacts)
        #timeseries.append([ts.frame, n_contacts])
timeseries= np.array(timeseries)

np.savetxt('dist_222_42_4a2j_NH2_OE2.txt',timeseries)
#print(timeseries.shape)
#sns.kdeplot(timeseries, bw_adjust=0.5)
#plt.title("Probability Distribution Function")
#plt.xlabel("Data values")
#plt.ylabel("Density")
#plt.show()


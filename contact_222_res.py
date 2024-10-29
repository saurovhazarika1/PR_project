import MDAnalysis as mda
from MDAnalysis.analysis import contacts

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob

import sys

trajectory_files = glob('2ovm/2ovm_strip_*.dcd')
print(trajectory_files)
topology_file = '2ovm/2ovm_strip.prmtop'


#sel_basic = "(resname ARG LYS) and (name NH1 NE NH2 NZ)"
#sel_acedic = "(resname ASP GLU) and (name OD1 OD2 OE1 OE2)"

sel_basic_1 = "protein and name CA"
sel_acedic_1 = "protein and name CA"

sel_basic_2 = "(resid 42) and (name NH1 NE NH2 NZ)"
sel_acedic_2 = "(resid 222) and (name OD1 OD2 OE1 OE2)"


#def contacts_within_cutoff(group_a, group_b, radius=8.0):
    
timeseries = []

for traj in trajectory_files:
    u = mda.Universe(topology_file, traj)
    for ts in u.trajectory:
        acedic = u.select_atoms(sel_basic_1)
        basic = u.select_atoms(sel_acedic_1)
        dist = np.array(contacts.distance_array(acedic.positions, basic.positions))
        #contact_res = np.sum(dist, axis=1)
        #is_symmetric = np.allclose(dist, dist.T)
        #print("Is the matrix symmetric?", is_symmetric)
        n_contacts = contacts.contact_matrix(dist, 8.0)
        n_contacts_222 = n_contacts[221,:]
        print(n_contacts_222[:5])
        timeseries.append(n_contacts_222)
        #timeseries = np.array(timeseries)
        #print(np.array(timeseries).shape)
        #timeseries.append([ts.frame, n_contacts])
#return np.array(timeseries)

timeseries = np.array(timeseries)
#timeseries_222 = np.sum(timeseries, axis=0)

np.savetxt('residue_222_contacts_t.dat', timeseries)

print(timeseries.shape)

#ca_1 = contacts_within_cutoff(sel_acedic_1, sel_basic_1, radius=8.0)

#print(ca_1.shape())
#ca_2 = contacts_within_cutoff(sel_acedic_2, sel_basic_2, radius=3.5)

#ca_1 = np.array(ca_1)
#ca_2 = np.array(ca_2)

#ca_H12_1 = ca_1[:,1]
#ca_H12_2 = ca_2[:,1]
#ca_H12_total = ca_H12_1 + ca_H12_2
#ca_H12_1 = ca_H12_1.reshape(-1,1)
#ca_H12_2 = ca_H12_2.reshape(-1,1)
#ca_H12_total = ca_H12_total.reshape(-1,1)

#ca_total = contacts_within_cutoff(sel_acedic, sel_basic, radius=3.5)

#ca_total = np.array(ca_total)
#ca_total = ca_total[:,1]
#ca_total = ca_total.reshape(-1,1)


#print(ca_1[:10])
#print(ca_2[:10])
#print(ca_H12_1[:10])
#print(ca_H12_2[:10])

#print(ca_1.shape)
#print(ca_2.shape)
#print(ca_H12_1.shape)
#print(ca_H12_2.shape)
#mean_val = np.mean(ca[:,1])
#print(mean_val)

#ca_1_df = pd.DataFrame(ca_1, columns=['Frame','# contacts'])
#ca_1_df.head()

#ca_2_df = pd.DataFrame(ca_2, columns=['Frame','# contacts'])
#ca_2_df.head()

#ca_df = pd.DataFrame(ca, columns=['# contacts'])
#ca_df.head()

#np.savetxt('saltbridge_217_41.dat',ca_H12_1)
#np.savetxt('saltbridge_222_42.dat',ca_H12_2)
#np.savetxt('saltbridge_total.dat',ca_H12_total)


#ca_df.plot(x='Frame')
#plt.ylabel('# salt bridges')
#plt.show()

#mean_val_2 = np.empty(len(ca[:,0]))
#mean_val_2.fill(mean_val)
#print(mean_val_2)

#plt.bar(ca[:,0], ca[:,1])
#plt.plot(ca[:,0], mean_val_2, color='r', label='mean')
#plt.xlabel('frames')
#plt.ylabel('number of salt bridges')
#plt.title('number of salt bridges formed by CDCA H12 residues')
#plt.legend()
#plt.show()

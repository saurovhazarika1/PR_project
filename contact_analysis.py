##This code is used to calculate the contacts between H12 region and L11-12 region##
##Written by Saurov Hazarika##

import mdtraj as md
import glob

#indir = '2ovm'
#top = f'{indir}/2ovm_strip.prmtop'
#trajs = glob(f'{indir}/2ovm_strip_*.dcd')

top = '2ovm/2ovm_strip.prmtop'
traj = '2ovm/2ovm_strip.dcd'

trajectory = md.load(traj, top=top)

#trajectory_list = []

#for traj_file in trajs:
#    for traj_chunk in md.iterload(traj_file, top=top):
#        trajectory_list.append(traj_chunk)


#trajectory = md.join(trajectory_list)

backbone_indices = trajectory.top.select('resid 10 to 215 and backbone')

print(backbone_indices)

#sidechain_indices = trajectory.top.select('resid 1 to 215 and sidechain')

H12_indices = trajectory.top.select('resid 226 to 246 and backbone')
#H12_indices_sidechain = trajectory.top.select('resid 226 to 246 and sidechain')

L11_12_indices = trajectory.top.select('resid 216 to 225 and backbone')
#L11_12_indices_sidechain = trajectory.top.select('resid 216 to 225 and sidechain')


print("Max backbone index:", max(backbone_indices))
print("Max H12 indices:", max(H12_indices))
print("Max L11_12 indices:", max(L11_12_indices))
print("Total number of atoms in trajectory:", trajectory.n_atoms)

import numpy as np

H12_pairs = [[i, j] for i in H12_indices for j in backbone_indices]

print(H12_pairs)
#H12_SC_pairs = np.array([[i, j] for i in H12_indices for j in sidechain_indices])

L11_12_pairs = [[i, j] for i in L11_12_indices for j in backbone_indices]
#L11_12_SC_pairs = np.array([[i, j] for i in L11_12_indices for j in sidechain_indices])

backbone_contact_distances_H12 = md.compute_contacts(trajectory, contacts=H12_pairs)
#sidechain_contact_distances_H12, _ = md.compute_contacts(trajectory, contacts=H12_SC_pairs)

backbone_contact_distances_L11 = md.compute_contacts(trajectory, contacts=L11_12_pairs)
#sidechain_contact_distances_L11, _ = md.compute_contacts(trajectory, contacts=L11_12_SC_pairs)

cutoff = 0.4  # nm

backbone_H12_contacts_within_cutoff = backbone_contact_distances_H12 < cutoff
#sidechain_H12_contacts_within_cutoff = sidechain_contact_distances_H12 < cutoff

backbone_L11_contacts_within_cutoff = backbone_contact_distances_L11 < cutoff
#sidechain_L11_contacts_within_cutoff = sidechain_contact_distances_L11 < cutoff

np.savetxt('backbone_H12_contacts.dat', backbone_H12_contacts_within_cutoff)
#np.savetxt('sidechain_H12_contacts.dat', sidechain_H12_contacts_within_cutoff)
np.savetxt('backbone_L11_contacts.dat', backbone_L11_contacts_within_cutoff)
#np.savetxt('sidechain_L11_contacts.dat', sidechain_L11_contacts_within_cutoff)


#num_backbone_H12_contacts_per_frame = np.sum(backbone_H12_contacts_within_cutoff, axis=1)
#num_sidechain_H12_contacts_per_frame = np.sum(sidechain_H12_contacts_within_cutoff, axis=1)

#num_backbone_L11_contacts_per_frame = np.sum(backbone_L11_contacts_within_cutoff, axis=1)
#num_backbone_L11_contacts_per_frame = np.sum(sidechain_L11_contacts_within_cutoff, axis=1)

#import matplotlib.pyplot as plt

#plt.plot(traj.time, num_backbone_contacts_per_frame, label='Backbone-Backbone Contacts')
#plt.plot(traj.time, num_sidechain_contacts_per_frame, label='Sidechain-Sidechain Contacts')
#plt.xlabel('Time (ps)')
#plt.ylabel('Number of Contacts')
#plt.legend()
#plt.show()



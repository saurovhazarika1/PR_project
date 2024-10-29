import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import seaborn as sns

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob

trajectory_files = glob('2ovm/2ovm_strip_*.dcd')
print(trajectory_files)
topology_file = '2ovm/2ovm_strip.prmtop'


protein_res = 'protein'
ligand_res = 'resname ASO'


def contacts_within_cutoff(group_a, group_b, radius=4.0):

    timeseries = []

    for traj in trajectory_files:
        u = mda.Universe(topology_file, traj)
        for ts in u.trajectory:
            protein = u.select_atoms(group_a)
            ligand = u.select_atoms(group_b)
            dist = contacts.distance_array(protein.positions, ligand.positions)
            n_contacts = contacts.contact_matrix(dist, radius).sum()
            timeseries.append([ts.frame, n_contacts])
    return np.array(timeseries)

ca_1 = contacts_within_cutoff(protein_res, ligand_res, radius=3.5)

ca_1 = np.array(ca_1)
ca_contacts = ca_1[:,1]

print(ca_contacts[0])
print(ca_contacts.shape)

np.savetxt('protein_ligand_dist.dat', ca_contacts)

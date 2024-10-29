##This code is used to do the PCA in MD trajectories using MDAnalysis for PR system##
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD
from MDAnalysis.analysis import pca, align

import warnings
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

top = sys.argv[1]
traj = sys.argv[2]

u = mda.Universe(top, traj)

aligner = align.AlignTraj(u, u, select='backbone',in_memory=True).run()

pc = pca.PCA(u, select='backbone',align=True, mean=None,n_components=2).run()

backbone = u.select_atoms('backbone')
n_bb = len(backbone)
#print('There are {} backbone atoms in the analysis'.format(n_bb))
#print(pc.p_components.shape)

print(pc.cumulated_variance[0])
print(pc.cumulated_variance[1])

transformed = pc.transform(backbone, n_components=3)
transformed.shape

df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(3)])
df['Time (ps)'] = df.index * u.trajectory.dt
df.head()

pc1 = pc.p_components[:, 0]
trans1 = transformed[:, 0]
projected = np.outer(trans1, pc1) + pc.mean.flatten()
coordinates = projected.reshape(len(trans1), -1, 3)

proj1 = mda.Merge(backbone)
proj1.load_new(coordinates, order="fac")

pdbtrj = 'pc1_newtraj_4a2j.pdb'

with mda.Writer(pdbtrj, multiframe=True, n_atoms=proj1.atoms.n_atoms) as PDB:
    for ts in proj1.trajectory:
        PDB.write(proj1.atoms)

pc2 = pc.p_components[:, 1]
trans2 = transformed[:, 1]
projected = np.outer(trans2, pc2) + pc.mean.flatten()
coordinates = projected.reshape(len(trans2), -1, 3)

proj2 = mda.Merge(backbone)
proj2.load_new(coordinates, order="fac")

pdbtrj = 'pc2_newtraj_4a2j.pdb'

with mda.Writer(pdbtrj, multiframe=True, n_atoms=proj1.atoms.n_atoms) as PDB:
    for ts in proj2.trajectory:
        PDB.write(proj2.atoms)








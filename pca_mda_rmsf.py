##This code calculates the RMSF along principal components##
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys

import MDAnalysis as mda
from MDAnalysis.tests.datafiles import PSF, DCD
from MDAnalysis.analysis import rms, pca, align

import warnings
# suppress some MDAnalysis warnings about writing PDB files
warnings.filterwarnings('ignore')

top = sys.argv[1]
traj = sys.argv[2]

u = mda.Universe(top, traj)

aligner = align.AlignTraj(u, u, select='backbone',in_memory=True).run()

pc = pca.PCA(u, select='backbone',align=True, mean=None,n_components=5).run()

backbone = u.select_atoms('backbone')
n_bb = len(backbone)
#print('There are {} backbone atoms in the analysis'.format(n_bb))
#print(pc.p_components.shape)

#print(pc.cumulated_variance[0])
#print(pc.cumulated_variance[1])
#print(pc.cumulated_variance[2])
#print(pc.cumulated_variance[3])
#print(pc.cumulated_variance[4])

transformed = pc.transform(backbone, n_components=3)
transformed.shape

#df = pd.DataFrame(transformed,
                  #columns=['PC{}'.format(i+1) for i in range(3)])
#df['Time (ps)'] = df.index * u.trajectory.dt
#df.head()

pc1 = pc.p_components[:, 0]
trans1 = transformed[:, 0]
projected = np.outer(trans1, pc1) + pc.mean.flatten()
coordinates = projected.reshape(len(trans1), -1, 3)

proj1 = mda.Merge(backbone)
proj1.load_new(coordinates, order="fac")

average = align.AverageStructure(proj1, proj1, select='protein and name CA', ref_frame=0).run()
ref = average.results.universe

aligner = align.AlignTraj(proj1, ref, select='protein and name CA', in_memory=True).run()

c_alphas_1 = proj1.select_atoms('protein and name CA')
#c_alphas_1_dbd = proj1.select_atoms('name CA and resid 1 to 75')
#c_alphas_1_hinge = proj1.select_atoms('name CA and resid 76 to 124')
#c_alphas_1_lbd = proj1.select_atoms('name CA and resid 125 to 356')

R1 = rms.RMSF(c_alphas_1).run()
#R1_dbd = rms.RMSF(c_alphas_1_dbd).run()
#R1_hinge = rms.RMSF(c_alphas_1_hinge).run()
#R1_lbd = rms.RMSF(c_alphas_1_lbd).run()


pc2 = pc.p_components[:, 1]
trans2 = transformed[:, 1]
projected = np.outer(trans2, pc2) + pc.mean.flatten()
coordinates = projected.reshape(len(trans2), -1, 3)

proj2 = mda.Merge(backbone)
proj2.load_new(coordinates, order="fac")

average = align.AverageStructure(proj2, proj2, select='protein and name CA', ref_frame=0).run()
ref = average.results.universe

aligner = align.AlignTraj(proj2, ref, select='protein and name CA', in_memory=True).run()

c_alphas_2 = proj2.select_atoms('protein and name CA')
#c_alphas_2_dbd = proj2.select_atoms('name CA and resid 1 to 75')
#c_alphas_2_hinge = proj2.select_atoms('name CA and resid 76 to 124')
#c_alphas_2_lbd = proj2.select_atoms('name CA and resid 125 to 356')

R2 = rms.RMSF(c_alphas_2).run()
#R2_dbd = rms.RMSF(c_alphas_2_dbd).run()
#R2_hinge = rms.RMSF(c_alphas_2_hinge).run()
#R2_lbd = rms.RMSF(c_alphas_2_lbd).run()

np.savetxt('pc1_4a2j_rmsf.dat',R1.results.rmsf)
#np.savetxt('pc1_m1_CDCA_dbd_rmsf.dat',R1_dbd.results.rmsf)
#np.savetxt('pc1_m1_CDCA_hinge_rmsf.dat',R1_dbd.results.rmsf)
#np.savetxt('pc1_m1_CDCA_lbd_rmsf.dat',R1_lbd.results.rmsf)

np.savetxt('pc2_4a2j_rmsf.dat',R2.results.rmsf)
#np.savetxt('pc2_m1_CDCA_dbd_rmsf.dat',R2_dbd.results.rmsf)
#np.savetxt('pc2_m1_CDCA_hinge_rmsf.dat',R2_hinge.results.rmsf)
#np.savetxt('pc2_m1_CDCA_lbd_rmsf.dat',R2_lbd.results.rmsf)


#print(len(R1_dbd.results.rmsf))
#print(len(R1_hinge.results.rmsf))
#print(len(R1_lbd.results.rmsf))
#print(len(R2_dbd.results.rmsf))
#print(len(R2_hinge.results.rmsf))
#print(len(R2_lbd.results.rmsf))


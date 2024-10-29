##This code calculates the pearson correlation of Sine and Cosine of dihedral angles with TICA1 and TICA2 componets##
##The dihedral angles can be calculated for the trajectories using a different code##
##-Saurov##

from __future__ import print_function
import pyemma
pyemma.__version__
import os
import matplotlib
import matplotlib.pyplot as plt
import pylab
matplotlib.rcParams.update({'font.size': 12})
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 15),
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'x-large',
         'ytick.labelsize':'x-large',
         'font.weight':'medium',
         'xtick.major.size':8,
         'ytick.major.size':8}
pylab.rcParams.update(params)
import pyemma.coordinates as coor
import pyemma.msm as msm
import pyemma.plots as pplt
import numpy as np
import pickle
import mdtraj as md
import itertools
from pyemma.coordinates import source
from pyemma.coordinates.data import CustomFeature
from numpy.linalg import norm
from IPython.display import display
import sys
from scipy import stats

a = sys.argv[2]
b = sys.argv[3]
c = sys.argv[4]
d = sys.argv[5]

phi_sin = np.genfromtxt(a)
phi_cos = np.genfromtxt(b)
psi_sin = np.genfromtxt(c)
psi_cos = np.genfromtxt(d)


dcd1 = sys.argv[1]
#psf1 = sys.argv[2]

#dcd2 = sys.argv[2]
#psf2 = sys.argv[4]

indir = '2ovm'
top =  indir+'/2ovm_strip.prmtop'
from glob import glob
trajs = glob ('2ovm/2ovm_strip_*.dcd')
print(trajs[38])

feat = pyemma.coordinates.featurizer(top)

traj1 = md.load(dcd1, top=top)

feat.add_minrmsd_to_ref(traj1)

region_1 = [i for i in range(22, 31)]
region_2 = [i for i in range(216, 226)]
pairs = []
for i in region_1:
    for j in region_2:
        pairs.append([i,j])
pairs = np.array(pairs)

feat.add_residue_mindist(pairs,periodic=False)

#feat = pyemma.coordinates.featurizer(psf2)

traj2 = md.load(dcd1, top=top)

feat.add_minrmsd_to_ref(traj2, ref_frame=-1)

reader = pyemma.coordinates.source(trajs, features=feat)

data_angle = reader.get_output()

tica_obj = coor.tica(data_angle, lag=100, dim=3, kinetic_map=False, var_cutoff=0.95, reversible=True)

tica_output = np.array(tica_obj.get_output())
print(tica_output.shape)
tica_concatenated = tica_output.flatten()
#tica_concatenated = tica_output.reshape(480000, 2)
print(tica_concatenated.shape)

#n_cluster = 1000
#clustering = coor.cluster_kmeans(tica_concatenated, k=n_cluster, max_iter=500)
#dtrajs = clustering.dtrajs

data = np.loadtxt('total_y_2ovm.out')
extended = data[:,1]
bridge = data[:,2]
three = data[:,3]
alpha = data[:,4]
pi = data[:,5]
turn = data[:,6]
bend = data[:,7]
helix = [0]*len(pi)
for i in range(len(pi)):
    b = three[i]
    c = alpha[i]
    d = pi[i]
    helix[i] = b+c+d
other = [0]*len(pi)
for i in range(len(pi)):
    b = bridge[i]
    c = turn[i]
    d = bend[i]
    other[i] = b+c+d
helix = np.array(helix)
beta = np.array(extended)
other = np.array(other)
x = np.vstack(tica_concatenated)
print(x.shape)
tica_concatenated = np.column_stack((x, helix))
tica_concatenated = np.column_stack((tica_concatenated, beta))
tica_concatenated = np.column_stack((tica_concatenated, other))
print(tica_concatenated.shape)


data_SB = np.loadtxt('saltbridge_217_41.dat')
tica_concatenated = np.column_stack((tica_concatenated, data_SB))

n_cluster = 1000
clustering = coor.cluster_kmeans(tica_concatenated, k=n_cluster, max_iter=500)
dtrajs = clustering.dtrajs

#t = np.vstack(tica_concatenated)[:,0]
itc1 = tica_concatenated[:,0]
itc2 = tica_concatenated[:,1]

np.savetxt('itc1.dat',itc1)
np.savetxt('itc2.dat',itc2)

corr_phi_sin_1 = []
corr_phi_cos_1 = []
corr_phi_sin_2 = []
corr_phi_cos_2 = []

corr_psi_sin_1 = []
corr_psi_cos_1 = []
corr_psi_sin_2 = []
corr_psi_cos_2 = []

for i in range(249):
    phi_sin_res = phi_sin[i,:]
    phi_cos_res = phi_cos[i,:]
    psi_sin_res = psi_sin[i,:]
    psi_cos_res = psi_cos[i,:]

    corr_phi_sin_1_res, p_value_1 = stats.pearsonr(phi_sin_res, itc1)
    corr_phi_cos_1_res, p_value_2 = stats.pearsonr(phi_cos_res, itc1)
    corr_phi_sin_2_res, p_value_3 = stats.pearsonr(phi_sin_res, itc2)
    corr_phi_cos_2_res, p_value_4 = stats.pearsonr(phi_cos_res, itc2)

    corr_psi_sin_1_res, p_value_5 = stats.pearsonr(psi_sin_res, itc1)
    corr_psi_cos_1_res, p_value_6 = stats.pearsonr(psi_cos_res, itc1)
    corr_psi_sin_2_res, p_value_7 = stats.pearsonr(psi_sin_res, itc2)
    corr_psi_cos_2_res, p_value_8 = stats.pearsonr(psi_cos_res, itc2)

    corr_phi_sin_1.append(corr_phi_sin_1_res)
    corr_phi_cos_1.append(corr_phi_cos_1_res)
    corr_phi_sin_2.append(corr_phi_sin_2_res)
    corr_phi_cos_2.append(corr_phi_cos_2_res)


    corr_psi_sin_1.append(corr_psi_sin_1_res)
    corr_psi_cos_1.append(corr_psi_cos_1_res)
    corr_psi_sin_2.append(corr_psi_sin_2_res)
    corr_psi_cos_2.append(corr_psi_cos_2_res)

np.savetxt('corr_phi_sin_1.txt', corr_phi_sin_1)
np.savetxt('corr_phi_cos_1.txt', corr_phi_cos_1)
np.savetxt('corr_phi_sin_2.txt', corr_phi_sin_2)
np.savetxt('corr_phi_cos_2.txt', corr_phi_cos_2)

np.savetxt('corr_psi_sin_1.txt', corr_phi_sin_1)
np.savetxt('corr_psi_cos_1.txt', corr_phi_cos_1)
np.savetxt('corr_psi_sin_2.txt', corr_phi_sin_2)
np.savetxt('corr_psi_cos_2.txt', corr_phi_cos_2)


#fig, ax = plt.subplots(1, 4, sharex=True, sharey=True)
#cc_x = clustering.clustercenters[:,0]
#cc_y = clustering.clustercenters[:,1]
#pyemma.plots.plot_free_energy(np.vstack(tica_concatenated)[:,0], np.vstack(tica_concatenated)[:,1], ax=ax[0], cbar=True)
#ax[1].scatter(cc_x,cc_y, c=clustering.clustercenters[:,3],cmap='gnuplot', marker='o', vmin=0,vmax=0.68)
#ax[1].set_xlabel('Helix')
#ax[2].scatter(cc_x, cc_y, c=clustering.clustercenters[:,4],cmap='gnuplot', marker='o', vmin=0,vmax=0.092)
#ax[2].set_xlabel('Beta')
#d = ax[3].scatter(cc_x, cc_y, c=clustering.clustercenters[:,5],cmap='gnuplot', marker='o', vmin=0,vmax=0.376)
#ax[3].set_xlabel('Other')
#d=ax[1].scatter(cc_x, cc_y, c=clustering.clustercenters[:,6],cmap='gnuplot', marker='o', vmin=0,vmax=4.0)
#ax[1].set_xlabel('SB')
#plt.colorbar(d, ax=ax[1])
#plt.tight_layout()
#plt.show()
#plt.savefig('ss_colored_sec_new.png',dpi=300)


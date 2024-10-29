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



indir = '2ovm'
top =  indir+'/2ovm_strip.prmtop'
from glob import glob
trajs = glob ('2ovm/2ovm_strip_*.dcd')
print(trajs)

traj = md.iterload(trajs, top=top)
feat = pyemma.coordinates.featurizer(top)
trajs = glob ('2ovm/2ovm_strip_*.dcd')
print(trajs)
reader = source(trajs, top=top)

A = feat.topology.select("resid 225 and name CA")

B = feat.topology.select("resid 241 and name CA")

C = feat.topology.select("resid 216 and name CA")

D = feat.topology.select("resid 181 and name CA")

#print(trajs)
def feature(traj: md.Trajectory):
    ang_list = []
    for i in range(len(traj)):
        vec1 = np.array(np.array(traj.xyz[i,B])-np.array(traj.xyz[i,A]))
        vec2 = np.array(np.array(traj.xyz[i,D])-np.array(traj.xyz[i,C]))
        angle = np.arccos(np.dot(vec1[0],vec2[0])/(norm(vec1[0])*norm(vec2[0])))
        angle = np.rad2deg(angle)
        ang_list.append(angle)
    ang_list = np.array(ang_list)
    ang_list_reshaped = ang_list.reshape(-1,1)
    ang_list_reshaped = ang_list_reshaped.astype('float32')
    return ang_list_reshaped
    print(ang_list_reshaped.shape)

#print(len(trajs))
#print(len(traj))
#print(feature(traj).shape)

my_feature = CustomFeature(feature, dim=1)

reader.featurizer.add_custom_feature(my_feature)

data_angle = reader.get_output()

tica_obj = coor.tica(data_angle, lag=100, kinetic_map=False, var_cutoff=0.95, reversible=True)

tica_output = tica_obj.get_output()
tica_concatenated = np.concatenate(tica_output)

n_cluster = 250
clustering = coor.cluster_kmeans(tica_output, k=n_cluster, max_iter=500)
dtrajs = clustering.dtrajs


#its = msm.timescales_msm(dtrajs, lags=150, nits=5, errors='bayes')
#pplt.plot_implied_timescales(its, units='steps', nits=10, linewidth=2)
#plt.title('implied timescale plot with errors')
#plt.xlabel('Lag Time (steps)')
#plt.ylabel('Timescale (steps)')
#plt.figure()
#plt.show()

M = msm.bayesian_markov_model(dtrajs, 500, reversible=True)

#plt.plot(M.timescales(), linewidth=0, marker='o')
#plt.xlabel('index')
#plt.ylabel('timescale')
#plt.show()

#pyemma.plots.plot_cktest(M.cktest(2))
#plt.show()

print('fraction of states used = ', M.active_state_fraction)
print('fraction of states used = ', M.active_state_fraction)

nstates = 4
M.pcca(nstates)
pcca_dist = M.metastable_distributions

indir = '2ovm'
top =  indir+'/2ovm_strip.prmtop'

pcca_samples = M.sample_by_distributions(M.metastable_distributions, 1)
#torsions_source = pyemma.coordinates.source(trajs, top=top, features=my_feature)
pyemma.coordinates.save_trajs(reader, pcca_samples, 
        outfiles=['pcca{}-samples.pdb'.format(n+1)
            for n in range(M.n_metastable)])

for i, s in enumerate(M.metastable_sets):
    print('pi_{} = {:f}'.format(i + 1, M.pi[s].sum()))

nstates = 4
mfpt = np.zeros((nstates, nstates))
for i in range(nstates):
    for j in range(nstates):
        mfpt[i, j] = M.mfpt(
                M.metastable_sets[i],
                M.metastable_sets[j])

from pandas import DataFrame
print('MFPT / ns: ')
display(DataFrame(np.round(mfpt, decimals=2), index=range(0, nstates), columns=range(0, nstates)).multiply(0.06))

A = M.metastable_sets[0]
B = np.concatenate(M.metastable_sets[1:])
print('MFPT 1 -> other: ({:6.1f} +- {:5.1f}) ns'.format(
    0.06*M.sample_mean('mfpt', A, B), 0.06*M.sample_std('mfpt', A, B)))

print('MFPT other -> 1: ({:.1f} +- {:5.1f}) ns'.format(
    0.06*M.sample_mean('mfpt', B, A), 0.06*M.sample_std('mfpt', B, A)))












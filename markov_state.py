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

#dcd1 = sys.argv[1]
dcd2 = sys.argv[1]

indir = '2ovm'
top =  indir+'/2ovm_strip.prmtop'
from glob import glob
trajs = glob ('2ovm/2ovm_strip_*.dcd')
print(trajs)


feat = pyemma.coordinates.featurizer(top)

traj1 = md.load(dcd2, top=top)

feat.add_minrmsd_to_ref(traj1, ref_frame=1)

region_1 = [i for i in range(22, 31)]
region_2 = [i for i in range(216, 226)]
pairs = []
for i in region_1:
    for j in region_2:
        pairs.append([i,j])
pairs = np.array(pairs)

feat.add_residue_mindist(pairs,periodic=False)

traj2 = md.load(dcd2, top=top)

feat.add_minrmsd_to_ref(traj2, ref_frame=1)

reader = pyemma.coordinates.source(trajs, features=feat) 

data_angle = reader.get_output()


tica_obj = coor.tica(data_angle, lag=100, dim=3, kinetic_map=False, var_cutoff=0.95, reversible=True)

#tica_output = np.array(tica_obj.get_output())
#print(tica_output.shape)
#tica_concatenated = tica_output.flatten()
#tica_concatenated = np.concatenate(tica_output)
#print(tica_concatenated.shape)

tica_output = tica_obj.get_output()
tica_concatenated = np.concatenate(tica_output)

#tica_concatenated = np.vstack(tica_concatenated)
#print(tica_concatenated.shape)

#tica_output = tica_obj.get_output()
#tica_concatenated = np.concatenate(tica_output)

n_cluster = 1000
clustering = coor.cluster_kmeans(tica_concatenated, k=n_cluster, max_iter=500)
dtrajs = clustering.dtrajs


its = msm.timescales_msm(dtrajs, lags=150, nits=5, errors='bayes')
pplt.plot_implied_timescales(its, units='steps', nits=10, linewidth=2)
plt.title('implied timescale plot with errors')
plt.xlabel('Lag Time (steps)')
plt.ylabel('Timescale (steps)')
#plt.figure()
plt.legend()
plt.savefig('its.png')
#plt.show()

M = msm.bayesian_markov_model(dtrajs, 500, reversible=True)

plt.plot(M.timescales(), linewidth=0, marker='o')
plt.xlabel('index')
plt.ylabel('timescale')
plt.legend()
plt.savefig('indextime.png')
#plt.show()

#print(M.timsecales())

plt.plot(M.timescales()[:-1]/M.timescales()[1:], linewidth=0, marker='o')
plt.plot(np.full(len(M.timescales()[:-1]/M.timescales()[1:]), 1.5))
#plt.title('Timescale separation')
#plt.xlabel('index')
#plt.ylabel('timescale separation')
#plt.figure()
plt.savefig('timescalesep.png')
#plt.show()

pyemma.plots.plot_cktest(M.cktest(2))
plt.show()

print('fraction of states used = ', M.active_state_fraction)
print('fraction of states used = ', M.active_state_fraction)

nstates = 2
M.pcca(nstates)
pcca_dist = M.metastable_distributions

indir = '2ovm'
top =  indir+'/2ovm_strip.prmtop'

print(M.n_metastable)

pcca_samples = M.sample_by_distributions(M.metastable_distributions, 1)
#torsions_source = pyemma.coordinates.source(trajs, top=top, features=my_feature)
#pyemma.coordinates.save_trajs(reader, pcca_samples, 
#        outfiles=['pcca{}-samples_1.pdb'.format(n+1)
#            for n in range(2)])

for i, s in enumerate(M.metastable_sets):
    print('pi_{} = {:f}'.format(i + 1, M.pi[s].sum()))

nstates = 2
mfpt = np.zeros((nstates, nstates))
for i in range(nstates):
    for j in range(nstates):
        mfpt[i, j] = M.mfpt(
                M.metastable_sets[i],
                M.metastable_sets[j])

from pandas import DataFrame
print('MFPT / ns: ')
display(DataFrame(np.round(mfpt, decimals=2), index=range(0, nstates), columns=range(0, nstates)).multiply(0.01))

mfpt = np.round(mfpt*0.01, decimals=2)


mfpt_new = mfpt*0.01
inverse_mfpt=np.zeros_like(mfpt_new)
nz=mfpt.nonzero()
inverse_mfpt[nz]=1.0/mfpt_new[nz]


A = M.metastable_sets[0]
B = M.metastable_sets[1]


print('MFPT 1 -> 2: ({:6.1f} +- {:5.1f}) ns'.format(0.01*M.sample_mean('mfpt', A, B), 0.01*M.sample_std('mfpt', A, B)))


print('MFPT 2 -> 1: ({:.1f} +- {:5.1f}) ns'.format(0.01*M.sample_mean('mfpt', B, A), 0.01*M.sample_std('mfpt', B, A)))



dtrajs_concatenated = np.concatenate(dtrajs)
M.pcca(2)
pcca_dist = M.metastable_distributions
metastable_traj = M.metastable_assignments[dtrajs_concatenated]
fig, ax = plt.subplots(figsize=(4, 4))
_,_, misc = pyemma.plots.plot_state_map(tica_concatenated[:, 0].T,tica_concatenated[:, 1].T, metastable_traj, ax=ax)
ax.set_xlabel('IC 1')
ax.set_ylabel('IC 2')
misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1) for i in range(2)])
fig.tight_layout()
plt.title('PCCA')
plt.savefig('PCCA.png', dpi=300)
#plt.show()


print(M.metastable_assignments.shape)
#tica=tica_output[0]
#fig, ax=plt.subplots(1,4,figsize=(20,10),sharex=True,sharey=True)
#ax[0].scatter(tica_concatenated[0:1000,0][M.metastable_assignments==0],tica_concatenated[0:1000,1][M.metastable_assignments==0])
#ax[1].scatter(tica_concatenated[0:1000,0][M.metastable_assignments==1],tica_concatenated[0:1000,1][M.metastable_assignments==1],c='red')
#ax[2].scatter(tica_concatenated[0:1000,0][M.metastable_assignments==2],tica_concatenated[0:1000,1][M.metastable_assignments==2],c='green')
#ax[3].scatter(tica_concatenated[:,0][M.metastable_assignments==3],tica_concatenated[:,1][M.metastable_assignments==3],c='purple')
#plt.show()

#nstate = 2
#state_labels=['1','2']
#highest_membership = M.metastable_distributions.argmax(1)
#coarse_state_centers = clustering.clustercenters[M.active_set[highest_membership]]
#fig, ax = plt.subplots(figsize=(10,10))
#metastable_traj = M.metastable_assignments[dtrajs_concatenated]
#_,_, misc = pyemma.plots.plot_state_map(tica_concatenated[:, 0].T,tica_concatenated[:, 1].T, metastable_traj, ax=ax)
#ax.set_xlabel('IC 1')
#ax.set_ylabel('IC 2')
#misc['cbar'].set_ticklabels([r'$\mathcal{S}_%d$' % (i + 1) for i in range(2)])

#pyemma.plots.plot_network(
#        inverse_mfpt,
#        pos=coarse_state_centers[:,:2],
#        arrow_label_format='%.1f ns',
#        arrow_labels=mfpt,
#        arrow_scale=1.5,
#        ax=ax,
#        state_labels=range(1, nstate+1),
#        size=10, alpha=1);
#ax.set_xlabel('tica1')
#ax.set_ylabel('tica2')
#ax.set_xlim(tica_concatenated[:,0].min(), tica_concatenated[:,0].max())
#ax.set_xlim(tica_concatenated[:,1].min(), tica_concatenated[:,1].max())
#plt.savefig("invMFPT.png",dpi=300)








fig,ax=plt.subplots(1,1,figsize=(10,10))
state_labels=['1','2']
highest_membership = M.metastable_distributions.argmax(1)
coarse_state_centers = np.array([[-2.25, -1.0], [0.03,-0.10]])
size = np.array([0.0008,0.0008])
colors=('b','r')
ax.set_xlabel('tica1')
ax.set_ylabel('tica2')
ax.set_xlim(tica_concatenated[:,0].min(), tica_concatenated[:,0].max())
ax.set_ylim(tica_concatenated[:,1].min(), tica_concatenated[:,1].max())
ax.scatter(tica_concatenated[0:1000,0][M.metastable_assignments==0],
tica_concatenated[0:1000,1][M.metastable_assignments==0],label='State 1',alpha=0.3)
ax.scatter(tica_concatenated[0:1000,0][M.metastable_assignments==1],
tica_concatenated[0:1000,1][M.metastable_assignments==1],c='red',label='State 2',alpha=0.3)
#ax.scatter(tica_concatenated[0:1000,0][M.metastable_assignments==2],
#tica_concatenated[0:1000,1][M.metastable_assignments==2],c='green',label='State 3',alpha=0.3)
#ax.scatter(tica[0:996,0][M.metastable_assignments==3],
#tica[0:996,1][M.metastable_assignments==3],c='purple',label='State 4',alpha=0.3)
pyemma.plots.plot_network(
    inverse_mfpt,
    pos=coarse_state_centers[:,:2],
    arrow_label_format='%.1f ns',
    arrow_labels=mfpt,
    arrow_scale=1.5,
    state_colors=colors,
    state_sizes=size,
    ax=ax,
    state_labels=range(1, nstates + 1),
    size=10,alpha=1);
plt.legend()
plt.savefig('transition_graph.png',dpi=300)
plt.show()







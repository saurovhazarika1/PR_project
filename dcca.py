import mdtraj as md

indir = 'cdca'
top =  indir+'/fxrfl_dry.pdb'
from glob import glob
trajs = glob('cdca/cdca_*.dcd')
print(trajs)

topology = md.load(top)

#trajectory_list = []

#for traj_file in trajs:
trajectory = md.iterload(trajs, top=top)
#trajectory_list.append(traj)

#trajectory = md.join(trajectory_list)

selection = topology.top.select('name CA')

cross_corr = md.compute_dynamical_cross_correlation(trajectory, atom_pairs=selection)

import matplotlib.pyplot as plt
import seaborn as sns

sns.heatmap(cross_corr, cmap='coolwarm', center=0)
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.title('Dynamic Cross-Correlation Matrix')
plt.show()



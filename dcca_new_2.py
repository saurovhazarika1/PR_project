import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

indir = 'cdca'
top = f'{indir}/test.pdb'
trajs = glob(f'{indir}/cdca_*.dcd')
print(trajs)
trajectory_list = []

for traj_file in trajs:
    for traj_chunk in md.iterload(traj_file, top=top):
        trajectory_list.append(traj_chunk)

trajectory = md.join(trajectory_list)

selection = trajectory.top.select('name CA')

coords = trajectory.atom_slice(selection).xyz

n_frames = coords.shape[0]
n_atoms = coords.shape[1]

mean_coords = np.mean(coords, axis=0)

cov_matrix = np.zeros((n_atoms, n_atoms))

for i in range(n_atoms):
    for j in range(i, n_atoms):
        cov_ij = np.mean(np.sum((coords[:, i, :] - mean_coords[i]) * (coords[:, j, :] - mean_coords[j]), axis=1))
        cov_matrix[i, j] = cov_ij
        cov_matrix[j, i] = cov_ij  # Since the covariance matrix is symmetric

std_devs = np.sqrt(np.diag(cov_matrix))
cross_corr_matrix = cov_matrix / np.outer(std_devs, std_devs)

sns.heatmap(cross_corr_matrix, cmap='coolwarm', center=0)
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.title('Dynamic Cross-Correlation Matrix')
plt.show()


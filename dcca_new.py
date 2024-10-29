import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from glob import glob

# Define input files and directory
indir = 'cdca'
top = f'{indir}/test.pdb'
trajs = glob(f'{indir}/cdca_strip.dcd')
print(trajs)

# Load the topology
topology = md.load(top)

# Initialize an empty list to store trajectory objects
trajectory_list = []

# Load each trajectory iteratively and append to the list
for traj_file in trajs:
    traj = md.iterload(traj_file, top=top)
    trajectory_list.append(traj)

# Concatenate the list of trajectory objects into a single trajectory
trajectory = md.join(trajectory_list)

# Select C-alpha atoms
selection = topology.top.select('name CA')

# Extract the coordinates of the selected atoms
coords = trajectory.atom_slice(selection).xyz

# Number of frames and number of selected atoms
n_frames = coords.shape[0]
n_atoms = coords.shape[1]

# Calculate the mean position of each atom over the trajectory
mean_coords = np.mean(coords, axis=0)

# Calculate the covariance matrix of the atomic displacements
cov_matrix = np.zeros((n_atoms, n_atoms))

for i in range(n_atoms):
    for j in range(i, n_atoms):
        cov_ij = np.mean(np.sum((coords[:, i, :] - mean_coords[i]) * (coords[:, j, :] - mean_coords[j]), axis=1))
        cov_matrix[i, j] = cov_ij
        cov_matrix[j, i] = cov_ij  # Since the covariance matrix is symmetric

# Calculate the dynamic cross-correlation matrix
std_devs = np.sqrt(np.diag(cov_matrix))
cross_corr_matrix = cov_matrix / np.outer(std_devs, std_devs)

# Plot the cross-correlation matrix
sns.heatmap(cross_corr_matrix, cmap='coolwarm', center=0)
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.title('Dynamic Cross-Correlation Matrix')
plt.show()


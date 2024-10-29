##This code produce contacts maps based on Linear Mutual Information##
##Written by Saurov Hazarika but some help were taken from ChatGPT##

import MDAnalysis as mda
import numpy as np

traj = 'traj/Dimer_strip.dcd'
top = 'traj/Dimer_strip.prmtop'


u = mda.Universe(top, traj)


ca_atoms = u.select_atoms('resid 1-357 and name CA')

num_atoms = len(ca_atoms)

avg_pos = np.zeros((num_atoms, 3))  
displacements = []

for ts in u.trajectory:
    positions = ca_atoms.positions
    avg_pos += positions
avg_pos /= len(u.trajectory)  

for ts in u.trajectory:
    positions = ca_atoms.positions
    displacement = positions - avg_pos  
    displacements.append(displacement)

displacements = np.array(displacements)

reshaped_displacements = displacements.reshape(len(displacements), -1)

covariance_matrix = np.cov(reshaped_displacements, rowvar=False)


def compute_lmi(covariance_matrix, num_atoms):
    lmi_matrix = np.zeros((num_atoms, num_atoms))

    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            variance_i = covariance_matrix[3*i:3*(i+1), 3*i:3*(i+1)].trace()  
            variance_j = covariance_matrix[3*j:3*(j+1), 3*j:3*(j+1)].trace()

            covariance_ij = covariance_matrix[3*i:3*(i+1), 3*j:3*(j+1)].trace()
            rho_ij = covariance_ij / np.sqrt(variance_i * variance_j)

            if np.abs(rho_ij) < 1:
                lmi_matrix[i, j] = -0.5 * np.log(1 - rho_ij**2)
                lmi_matrix[j, i] = lmi_matrix[i, j]  

    return lmi_matrix

lmi_matrix = compute_lmi(covariance_matrix, num_atoms)


import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

#cmap = plt.get_cmap('coolwarm', 256)
plt.figure(figsize=(10, 8))
#sns.heatmap(lmi_matrix, cmap='coolwarm', square=True, cbar=True)
sns.heatmap(lmi_matrix, cmap='coolwarm', square=True, cbar=True, cbar_kws={'ticks': np.linspace(0, 3, 20)})
#plt.colorbar(ticks=np.linspace(0, 1, 10))

#custom_cmap = ListedColormap(colors)

# Plot the heatmap with the custom colormap
#sns.heatmap(lmi_matrix, cmap=custom_cmap, square=True, cbar=True)

plt.title('LMI-based Cross-Correlation Matrix')
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.show()

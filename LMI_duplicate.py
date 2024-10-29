import MDAnalysis as mda
import numpy as np

traj = '4a2j/4a2j_strip.dcd'
top = '4a2j/4a2j_strip.prmtop'


# Load your trajectory (e.g., PDB and DCD/XT files)
u = mda.Universe(top, traj)


# Select the C-alpha atoms
ca_atoms = u.select_atoms('name CA')

# Get the number of selected atoms (for building covariance matrices)
num_atoms = len(ca_atoms)


# Initialize the array to store average positions and displacements
avg_pos = np.zeros((num_atoms, 3))  # Average position for each atom
displacements = []

# Loop over trajectory frames to compute the average positions
for ts in u.trajectory:
    positions = ca_atoms.positions
    avg_pos += positions
avg_pos /= len(u.trajectory)  # Compute the average position

# Now, calculate the displacement vectors for each frame
for ts in u.trajectory:
    positions = ca_atoms.positions
    displacement = positions - avg_pos  # Displacement from average
    displacements.append(displacement)

displacements = np.array(displacements)


# Reshape the displacements for covariance calculation
reshaped_displacements = displacements.reshape(len(displacements), -1)

# Compute the covariance matrix (atoms x atoms)
covariance_matrix = np.cov(reshaped_displacements, rowvar=False)


def compute_lmi(covariance_matrix, num_atoms):
    # Initialize an LMI matrix (for simplicity, consider only inter-atomic LMI, not xyz components separately)
    lmi_matrix = np.zeros((num_atoms, num_atoms))

    # Loop through all pairs of atoms to calculate LMI
    for i in range(num_atoms):
        for j in range(i+1, num_atoms):
            # Variance of atom i and j
            variance_i = covariance_matrix[3*i:3*(i+1), 3*i:3*(i+1)].trace()  # Sum of variances along x, y, z
            variance_j = covariance_matrix[3*j:3*(j+1), 3*j:3*(j+1)].trace()

            # Covariance between i and j (sum of covariances along x, y, z)
            covariance_ij = covariance_matrix[3*i:3*(i+1), 3*j:3*(j+1)].trace()

            # Correlation coefficient
            rho_ij = covariance_ij / np.sqrt(variance_i * variance_j)

            # Compute LMI
            if np.abs(rho_ij) < 1:  # To avoid log(0) errors
                lmi_matrix[i, j] = -0.5 * np.log(1 - rho_ij**2)
                lmi_matrix[j, i] = lmi_matrix[i, j]  # Symmetric matrix

    return lmi_matrix

lmi_matrix = compute_lmi(covariance_matrix, num_atoms)
print(lmi_matrix.shape)

np.savetxt('lmi_matrix_4a2j.dat',lmi_matrix)

import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import ListedColormap

#cmap = plt.get_cmap('coolwarm', 256)
# Plot the LMI matrix as a heatmap
plt.figure(figsize=(10, 8))
#sns.heatmap(lmi_matrix, cmap='coolwarm', square=True, cbar=True)
sns.heatmap(lmi_matrix, cmap='coolwarm', square=True, cbar=True)
#plt.colorbar(ticks=np.linspace(0, 1, 10))

#custom_cmap = ListedColormap(colors)

# Plot the heatmap with the custom colormap
#sns.heatmap(lmi_matrix, cmap=custom_cmap, square=True, cbar=True)

#plt.title('LMI-based Cross-Correlation Matrix')
#plt.xlabel('Residue Index')
#plt.ylabel('Residue Index')
plt.show()

import numpy as np
import sys

a = sys.argv[1]
b = sys.argv[2]

data_mono = np.genfromtxt(a)
data_dim = np.genfromtxt(b)

lmi_diff_matrix = data_mono - data_dim

import matplotlib.pyplot as plt
import seaborn as sns

plt.figure(figsize=(10, 8))
#sns.heatmap(lmi_diff_matrix, cmap='coolwarm', center=0, square=True)
sns.heatmap(lmi_diff_matrix, cmap='coolwarm', square=True, cbar=True)
#plt.title('Difference in LMI-based Cross-Correlation Matrix')
#plt.xlabel('Residue Index')
#plt.ylabel('Residue Index')
plt.show()

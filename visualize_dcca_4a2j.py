import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

dccm = np.loadtxt('dcm_4a2j.dat')

sns.heatmap(dccm, cmap='coolwarm', center=0)
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.title('Dynamic Cross-Correlation Matrix 4a2j')
#plt.show()
plt.savefig('DCCM_4a2j.png')

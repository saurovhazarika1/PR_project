import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

dccm = np.loadtxt('dcm.dat')

sns.heatmap(dccm, cmap='coolwarm', center=0)
plt.xlabel('Residue Index')
plt.ylabel('Residue Index')
plt.title('Dynamic Cross-Correlation Matrix')
#plt.show()
plt.savefig('DCCM_2ovm.png')

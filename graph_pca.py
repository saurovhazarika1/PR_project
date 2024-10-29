import numpy as np
import sys
import matplotlib.pyplot as plt

a = sys.argv[1]
b = sys.argv[2]

data1 = np.genfromtxt(a)
data2 = np.genfromtxt(b)

plt.plot(data1, label='PC1')
plt.plot(data2, label='PC2')
plt.legend()
plt.savefig('PC_RMSF_4a2j.png')
plt.show()

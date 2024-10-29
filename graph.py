import numpy as np
import sys
import matplotlib.pyplot as plt

a = sys.argv[1]
data = np.genfromtxt(a)
data_mean = np.mean(data)

#data_mean_arr = np.ones(len(data))
#data_mean_arr = data_mean_arr*data_mean

plt.plot(data)
#plt.plot(data_mean_arr, label='mean')
#plt.legend()
plt.savefig('saltbridge_222_42_22.png')
plt.show()

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import seaborn as sns

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from glob import glob

import sys

a=sys.argv[1]
b=sys.argv[2]

data_1 = np.genfromtxt(a)
data_2 = np.genfromtxt(b)

sns.kdeplot(data_1, bw_adjust=0.5, label='Inactive')
sns.kdeplot(data_2, bw_adjust=0.5, label='active')
plt.legend()
plt.show()


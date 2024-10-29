import numpy as np


data = np.loadtxt('total_y_2ovm.out')
extended = data[:,1]
bridge = data[:,2]
three = data[:,3]
alpha = data[:,4]
pi = data[:,5]
turn = data[:,6]
bend = data[:,7]
helix = [0]*len(pi)
for i in range(len(pi)):
    b = three[i]
    c = alpha[i]
    d = pi[i]
    helix[i] = b+c+d
other = [0]*len(pi)
for i in range(len(pi)):
    b = bridge[i]
    c = turn[i]
    d = bend[i]
    other[i] = b+c+d
helix = np.array(helix)
beta = np.array(extended)
other = np.array(other)

print(max(helix))
print(max(beta))
print(max(other))

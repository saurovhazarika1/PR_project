import numpy as np
import sys

a = sys.argv[1]
b = sys.argv[2]
c = sys.argv[3]
d = sys.argv[4]


phi_sin_1 =np.absolute(np.genfromtxt(a)[1:])
phi_sin_2 =np.absolute(np.genfromtxt(b)[1:])
phi_cos_1 =np.absolute(np.genfromtxt(c)[1:])
phi_cos_2 =np.absolute(np.genfromtxt(d)[1:])

CS_score = phi_sin_1 + phi_sin_2 + phi_cos_1 + phi_cos_2

print(len(CS_score))

max_value = np.max(CS_score)
index = np.argmax(CS_score)

print('max_value = ', max_value)
print('index = ', index)

np.savetxt('CS_score_psi.txt', CS_score)

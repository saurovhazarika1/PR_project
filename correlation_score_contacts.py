import numpy as np
import sys

a = sys.argv[1]
b = sys.argv[2]


contacts_1 =np.absolute(np.genfromtxt(a)[1:])
contacts_2 =np.absolute(np.genfromtxt(b)[1:])

CS_score = contacts_1 + contacts_2 

print(len(CS_score))

max_value = np.nanmax(CS_score)
index = np.nanargmax(CS_score)

print('max_value = ', max_value)
print('index = ', index)

np.savetxt('CS_score_contacts_222.txt', CS_score)

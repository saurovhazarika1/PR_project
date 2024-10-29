import numpy as np
import sys
import matplotlib.pyplot as plt

a = sys.argv[1]
b = sys.argv[2]

data_217_41 = np.genfromtxt(a)
data_222_42 = np.genfromtxt(b)

matrix = np.zeros((5,5))

for i in range(len(data_217_41)):
    if data_217_41[i]==0:
        #pass
        if data_222_42[i]==0:
        #    pass
            matrix[0,0] += 1
        elif data_222_42[i]==1:
            matrix[0,1] += 1
        elif data_222_42[i]==2:
            matrix[0,2] += 1
        elif data_222_42[i]==3:
            matrix[0,3] += 1
        else:
            matrix[0,4] += 1

    elif data_217_41[i]==1:
        if data_222_42[i]==0:
            #pass
            matrix[1,0] += 1
        elif data_222_42[i]==1:
            matrix[1,1] += 1
        elif data_222_42[i]==2:
            matrix[1,2] += 1
        elif data_222_42[i]==3:
            matrix[1,3] += 1
        else:
            matrix[1,4] += 1

    elif data_217_41[i]==2:
        if data_222_42[i]==0:
            #pass
            matrix[2,0] += 1
        elif data_222_42[i]==1:
            matrix[2,1] += 1
        elif data_222_42[i]==2:
            matrix[2,2] += 1
        elif data_222_42[i]==3:
            matrix[2,3] += 1
        else:
            matrix[2,4] += 1

    elif data_217_41[i]==3:
        if data_222_42[i]==0:
            #pass
            matrix[3,0] += 1
        elif data_222_42[i]==1:
            matrix[3,1] += 1
        elif data_222_42[i]==2:
            matrix[3,2] += 1
        elif data_222_42[i]==3:
            matrix[3,3] += 1
        else:
            matrix[3,4] += 1

    else:
        if data_222_42[i]==0:
            #pass
            matrix[4,0] += 1
        elif data_222_42[i]==1:
            matrix[4,1] += 1
        elif data_222_42[i]==2:
            matrix[4,2] += 1
        elif data_222_42[i]==3:
            matrix[4,3] += 1
        else:
            matrix[4,4] += 1

plt.imshow(matrix, cmap='viridis', interpolation='none')  # 'viridis' is a color map
plt.colorbar()  # Add a color bar to show the scale
#plt.title("Heatmap of a 5x5 Matrix")
plt.show()

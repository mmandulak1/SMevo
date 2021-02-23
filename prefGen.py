import sys 
import random
import numpy as np

class Pairing:
    def __init__(self,x=-1,y=-1):
        self.pair = [x,y]

N = int(sys.argv[1])
matrix = np.empty((N,N), dtype=Pairing)
for i in range(0,N):
    for j in range(0,N):
        matrix[i][j]=Pairing()
outfile = open('prefList.txt','w')

for i in range(0,N):
        value = 1
        while value != N+1:
            coord = random.randint(1,N)-1
            if matrix[i][coord].pair[0] is -1:
                matrix[i][coord].pair[0] = value
                value += 1
        value = 1
        while value != N+1:
            coord = random.randint(1,N)-1
            if matrix[coord][i].pair[1] is -1:
                matrix[coord][i].pair[1] = value
                value += 1
for i in range(0,N):
    for j in range(0,N):
        if j == N-1:
            x = str(matrix[i][j].pair[0]) + "," + str(matrix[i][j].pair[1])
        else:
            x = str(matrix[i][j].pair[0]) + "," + str(matrix[i][j].pair[1]) + " "
        outfile.write(x)
    outfile.write('\n')
outfile.close()

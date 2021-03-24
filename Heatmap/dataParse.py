import matplotlib.pyplot as plt
import numpy as np
import sys


N = int(sys.argv[1])
matrix = np.zeros((N,N))
print(sys.argv[1])
data = open("slurm-3771.out", "r")
i=0
entry = ""
while True:
    x = data.read(1)
    if i < 5:
        i+=1
        continue
    if not x or x == "":
        break
    if x == '\n':
        continue
    if x == ',':
        if entry == "":
            continue
        index1 = int(entry)
        #print(entry)
        entry = ""
    elif x == " ":
        if entry == "":
            continue
        index2 = int(entry)
        #print(entry)
        matrix[index1][index2] += 1
        index1 = -1
        index2 = -1
        entry = ""
    else:
        entry += x
        
for i in range(N):
    for j in range(N):
        print(str(matrix[i][j])," ")



data.close()

import numpy as np 
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])
preflist = []
for i in range(0,N):
    new = []
    new2 = []
    for j in range(0,N):
        for k in range(0,2):
            new2.append(0)
        new.append(new2)
    preflist.append(new)

preflisttxt = open("PrefLists.txt.1.11.30","r")
entry = ""
i=0
j=0
while True:
    x = preflisttxt.read(1)
    #print(x)
    if not x or i==N:
        break
    if x == ',':
        if entry == "":
            continue
        leftval = int(entry)
        entry = ""
    elif x == " ":
        if entry == "":
            continue
        rightval = int(entry)
        entry = ""
        preflist[i][j]=[leftval,rightval]
        j+=1
    elif x == '\n':
        if entry == "":
            continue
        rightval = int(entry)
        entry = ""
        preflist[i][j]=[leftval,rightval]
        i+=1
        j=0
    else:
        entry += x
    if j==N:
        i+=1
        j=0



matrix = np.zeros((N,N))
heatvals = open("heatmapvals.txt","r")
flag = 0
i=0
j=0

while True:
    x = heatvals.readline()
    if flag == 0:
        flag+=1
        continue
    matrix[i][j] = float(x)
    j+=1
    if j==N:
        i+=1
        j=0
    if i==N:
        break

totalSum = 0
for i in range(N):
    for j in range(N):
        totalSum += matrix[i][j]

for i in range(N):
    for j in range(N):    
        matrix[i][j] = (matrix[i][j]/totalSum)*100

fig, ax = plt.subplots()
im = ax.imshow(matrix)

cbar = ax.figure.colorbar(im,ax=ax)
cbar.ax.set_ylabel("Percent Inclusion",rotation=-90,va="bottom")


ax.tick_params(top=True,bottom=False,labeltop=True,labelbottom=False)
for i in range(0,N):
    for j in range(0,N):
        title = str(preflist[i][j][0]) + "," + str(preflist[i][j][1])
        text = ax.text(j, i,title, ha="center", va="center", color="w")
ax.set_title("Pairing Inclusion for Random Cyclic Sample of N=11")
fig.tight_layout()
plt.show()
heatvals.close()
preflisttxt.close()
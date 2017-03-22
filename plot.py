import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
import time

npivot = 3
f = open("Results.txt", 'r')

p = []

i = 1
w = []
ans = []
for line in f:
    if ',' in line:
        y = [x for x in line[:-1].split(',')]
        z = [float(i) for i in y]
        z = tuple(z)
        w.append(z)
    else:
        l = line[:-1]
        l = float(l)
        w.append(l)
    if(i%(2*(npivot+2)+3) == 0):
        ans.append(w)
        w = []
    i+=1

numTimeSteps = int(ans[-1][0])+1
timestep = []
for i in range(numTimeSteps):
    timestep.append([])
for el in ans:
    timestep[int(el[0])].append(el)




def animate(timestep, i):
    plt.ion
    fig, ax = plt.subplots()
    npivot = 3
    npoint = 5
    currentNumPart = int(timestep[i][-1][1])
    for k in range(currentNumPart + 1):
        xc = []
        yc = []
        for j in range(npoint):
            xpoint = timestep[i][k][3+j][0]
            ypoint = timestep[i][k][3+j][1]
            xc.append(xpoint)
            yc.append(ypoint)
            circle1 = plt.Circle((xpoint, ypoint), timestep[i][k][2]/2, fill=False)
            ax.add_artist(circle1)
        plt.plot(xc, yc, '*-')
    plt.xlim([-10, 50])
    plt.ylim([-20, 20])
    plt.show()
    return None

for i in range(len(timestep)):
    if(i % 100 == 0):
        animate(timestep, i)

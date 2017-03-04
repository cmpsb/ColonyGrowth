import numpy as np
import matplotlib.pyplot as plt

inp = [1, (0, 0), (1, 2.5), (2, 4), (3, 2.5), (4, 0)]

npivot = 3

f = open("Results.txt", 'r')

p = []
numParticles = 3

for particle in range(numParticles):
    q = []
    for i in range(npivot+3):
        if i == 0:
            ans = float(f.readline())
        else:
            line = f.readline()
            line_list = line[:-1].split(",")
            ans = (line_list[0], line_list[1])
        q.append(ans)
    p.append(q)

print(p)

fig, ax = plt.subplots()
ax.set_xlim((-10, 10))
ax.set_ylim((-10, 10))
colors = ['r', 'g', 'b']
for i in range(numParticles):
    inp = p[i]
    coords = inp[1:]
    x = []
    y = []
    for j, c in enumerate(coords):
        x.append(c[0])
        y.append(c[1])
        circle1 = plt.Circle((x[j], y[j]), inp[0]/2, color=colors[i], fill=False)
        ax.add_artist(circle1)
    plt.plot(x, y, colors[i]+'x-')
ax.legend(['Daughter Particle 1', 'Daughter Particle 2', 'Mother Particle'])
plt.show()
del(fig)
del(ax)
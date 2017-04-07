import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

##========================================
npivot = 3
##========================================

fig, ax = plt.subplots()
x = np.arange(0, 2*np.pi, 0.01)
numparticles = 5
myLines = []
for j in range(numparticles):
    line, = ax.plot(0, 0)
    myLines.append(line)
for k in range((npivot+2)*numparticles):
    circle = patches.Circle((0, 0), 0, color='b', fill=False)
    ax.add_artist(circle)
    myLines.append(circle)

class ParticleData(object):
    def __init__(self):
        self.ts = 0
        self.ID = 0
        self.D = 0
        self.positions = []
        self.forces = []

    def setPosition(self, inp):
        self.positions.append(inp)

    def setForce(self, inp):
        self.forces.append(inp)
    def __str__(self):
        stringy = ["Time step is " + str(self.ts) + "\n"]
        stringy.append("Particle ID is " + str(self.ID) + "\n")
        stringy.append("Diameter is " + str(self.D) + "\n")
        for i, pos in enumerate(self.positions):
            stringy.append("Position {} has coordinates ({},{}) \n".format(i, self.positions[i][0], self.positions[i][1]))
        ans = "".join(stringy)
        return ans

def formatData(data_in):
    '''Make data per time step per particle'''
    data = []
    for line in data_in:
        timestep = []
        particles = line.split(";")
        for part in particles[:-1]:
            properties = part.split(" ")
            properties = properties[:-1]
            myPart = ParticleData()
            for i, prop in enumerate(properties):
                if i < 1:
                    myPart.ts = int(prop)
                elif i < 2:
                    myPart.ID = int(prop)
                elif i < 3:
                    myPart.D = float(prop)
                elif i < npivot + 5:
                    floats = [float(i) for i in prop.split(',')]
                    floats = tuple(floats)
                    myPart.setPosition(floats)
                elif i < 2*(npivot + 2) + 3:
                    floats = [float(i) for i in prop.split(',')]
                    floats = tuple(floats)
                    myPart.setForce(floats)
            timestep.append(myPart)
        data.append(timestep)
    return data

def animate(i):
    for j, particle in enumerate(data[i]):
        xs = []
        ys = []
        for k, point in enumerate(particle.positions):
            xs.append(point[0])
            ys.append(point[1])
            foo = (point[0], point[1])
            myLines[numparticles + j * (npivot+2) + k].center = foo
            myLines[numparticles + j * (npivot+2) + k].radius = particle.D/2
            myLines[numparticles + j * (npivot + 2) + k].set_color(myLines[j].get_color())
        myLines[j].set_xdata(xs)
        myLines[j].set_ydata(ys)
    return line,




f = open("Results.txt", 'r')
global data
data = formatData(f)
ax.set_xlim([-10, 10])
ax.set_ylim([-10, 10])
ani = animation.FuncAnimation(fig, animate, np.arange(19700,20420), interval=10)
ani.save(filename='myvid.mp4')	


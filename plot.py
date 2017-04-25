import matplotlib as mpl
import matplotlib.lines as lines
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import time

##========================================
npivot = 3
numparticles = 200
D = 1.0
range_x = [-10, 10]
range_y = [-10, 10]
##========================================
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


filename = input("Enter file name: ")
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
myLines = []
for j in range(numparticles):
    line, = ax.plot(0, 0)
    line.set_marker("o")
    line.set_markersize(0)
    line.set_solid_capstyle('round')
    myLines.append(line)

myShapes = []
# Stored in the form [rectangles first, then circles]
for j in range(numparticles):
    particleShapes = []
    for n in range(npivot+1):
        rect = patches.Rectangle((0, 0), width=0, height=0, angle=0.0)
        ax.add_artist(rect)
        particleShapes.append(rect)
    for n in range(npivot+2):
        circ = patches.Circle((0, 0), radius=0)
        ax.add_artist(circ)
        particleShapes.append(circ)
    myShapes.append(particleShapes)

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

def chooseColor():
    return (0, 0, 1)

def linewidth_from_data_units(linewidth, axis, reference='y'):
    fig = axis.get_figure()
    if reference == 'x':
        length = fig.bbox_inches.width * axis.get_position().width
        value_range = np.diff(axis.get_xlim())[0]
    elif reference == 'y':
        length = fig.bbox_inches.width * axis.get_position().width
        value_range = np.diff(axis.get_ylim())[0]
    length *= 72
    return linewidth * (length / value_range)

def makeRectangle(c1, c2, D):
    h = D
    v = np.array([(c2[0] - c1[0]), (c2[1] - c1[1])])
    w = np.linalg.norm(v)
    r = np.array([v[1], -v[0]])
    r *= D/(2*np.linalg.norm(r))
    corner = (c1[0] + r[0], c1[1] + r[1])
    theta = np.arctan2(v[1], v[0])
    return (corner, w, h, theta)

def animate(i):
    for j, particle in enumerate(data[i]):
        print(j)
        xs = []
        ys = []
        for point in particle.positions:
            xs.append(point[0])
            ys.append(point[1])
        color = chooseColor()
        myLines[j].set_color(color)
        myLines[j].set_xdata(xs)
        myLines[j].set_ydata(ys)
        myLines[j].set_linewidth(linewidth_from_data_units(1, ax, reference='y'))
        myLines[j].set_markersize(linewidth_from_data_units(0.5, ax, reference='y'))
    return line,

#takes 10E-4 sec, but still very slow...?
def neoAnimate(i):
    for j, particle in enumerate(data[i]):
        #set data rectangles
        for k in range(npivot+1):
            rectData = makeRectangle(particle.positions[k], particle.positions[k+1], D)
            myShapes[j][k].set_xy(rectData[0])
            myShapes[j][k].set_width(rectData[1])
            myShapes[j][k].set_height(rectData[2])
            transf = mpl.transforms.Affine2D().rotate_around(rectData[0][0], rectData[0][1], rectData[3]) + ax.transData
            myShapes[j][k].set_transform(transf)
        for l in range(npivot+1, 2*npivot + 3):
            myShapes[j][l].set_radius(D / 2)
            myShapes[j][l].center = (particle.positions[l - npivot - 1][0], particle.positions[l - npivot - 1][1])
    t1 = time.time()
    print(t1)
    return myShapes

f = open(filename, 'r')
global data

data = formatData(f)

ax.set_xlim(range_x)
ax.set_ylim(range_y)
ani = animation.FuncAnimation(fig, neoAnimate, np.arange(0, len(data)), interval=200, repeat=False)
plt.show()

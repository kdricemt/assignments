%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
from matplotlib.mlab import griddata

dir_trajectory = "./datafile/trajectory/"
nline = 170
yscale = 0.2

plt.rcParams["font.size"] = 30
fig = plt.figure(figsize=(60,20))

def plot_line():
    global nline,yscale
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.xlim(-10,15)
    plt.ylim(-10*yscale-2.0,10*yscale+2.0)
    plt.plot([-0.5,-0.5],[-0.5,0.5],color="black",linewidth=2.5)
    plt.plot([-0.5,0.5],[0.5,0.5],color="black",linewidth=2.5)
    plt.plot([0.5,0.5],[0.5,-0.5],color="black",linewidth=2.5)
    plt.plot([0.5,-0.5],[-0.5,-0.5],color="black",linewidth=2.5)


    for line in range(0,nline):
        filename = dir_trajectory + "trajectory" + str(line) + ".txt"
        data = np.loadtxt(filename)
        x = data[:,0]
        y = data[:,1]
        plt.plot(x,y,color = "black",linestyle="dashed")
        print("line:",line)

    filename_save = "./fig/trajectory.png"
    plt.savefig(filename_save)
    plt.show()

def plot_points():
    global nline,yscale
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.xlim(-10,15)
    plt.ylim(-10*yscale-2.0,10*yscale+2.0)
    plt.plot([-0.5,-0.5],[-0.5,0.5],color="black",linewidth=2.5)
    plt.plot([-0.5,0.5],[0.5,0.5],color="black",linewidth=2.5)
    plt.plot([0.5,0.5],[0.5,-0.5],color="black",linewidth=2.5)
    plt.plot([0.5,-0.5],[-0.5,-0.5],color="black",linewidth=2.5)

    for line in range(0,nline):
        filename = dir_trajectory + "trajectory" + str(line) + ".txt"
        data = np.loadtxt(filename)
        x = data[:,0]
        y = data[:,1]
        x2 = x[0:nline:10]
        y2 = y[0:nline:10]
        plt.scatter(x2,y2,marker='o',color='black')
        print("line:",line)

    filename_save = "./fig/tracer_points.png"
    plt.savefig(filename_save)
    plt.show()

if __name__ == "__main__":
    #plot_line()
    plot_points()

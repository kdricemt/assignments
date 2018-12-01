%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
from matplotlib.mlab import griddata

nline = 150

#parameters for animation
dstep = 200
nstep = 5000
nframe = int((nstep-10)/dstep) - 1;
yscale = 0.2


plt.rcParams["font.size"] = 30
fig = plt.figure(figsize=(60,20))
dir_streamline = "./datafile/streamline/"

def init():
    print("calculating step:",10)
    for line in range(0,nline):
        step = 10
        filename = dir_streamline + str(step) + "_" + str(line) + ".txt"
        data = np.loadtxt(filename)
        x = data[:,0]
        y = data[:,1]

        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xlim(-10,20)
        plt.ylim(-10*yscale-2.0,10*yscale+2.0)
        plt.plot(x,y,color="black",linestyle="dashed")
        plt.plot([-0.5,-0.5],[-0.5,0.5],color="black",linewidth=2.5)
        plt.plot([-0.5,0.5],[0.5,0.5],color="black",linewidth=2.5)
        plt.plot([0.5,0.5],[0.5,-0.5],color="black",linewidth=2.5)
        plt.plot([0.5,-0.5],[-0.5,-0.5],color="black",linewidth=2.5)

        txt = 'Steps:' + str(10) + ', ' + 'Time:' + str(0*100/50)
        plt.text(20,7.5,txt)


def animate(n):
    plt.cla() #delete previent data
    step = (n + 1)*dstep + 10
    print("calculating step:",step)

    plt.xlabel('X')
    plt.ylabel('Y')
    plt.xlim(-10,20)
    plt.ylim(-10*yscale-2.0,10*yscale+2.0)
    plt.plot([-0.5,-0.5],[-0.5,0.5],color="black",linewidth=2.5)
    plt.plot([-0.5,0.5],[0.5,0.5],color="black",linewidth=2.5)
    plt.plot([0.5,0.5],[0.5,-0.5],color="black",linewidth=2.5)
    plt.plot([0.5,-0.5],[-0.5,-0.5],color="black",linewidth=2.5)


    for line in range(0,nline):
        filename = dir_streamline + str(step) + "_" + str(line) + ".txt"
        data = np.loadtxt(filename)
        x = data[:,0]
        y = data[:,1]
        plt.plot(x,y,color = "black",linestyle="dashed")
        #print("step,line:",step," ",line,"(x,y): ",x,y,"   ",filename
    txt = 'Steps:' + str(step) + ', ' + 'Time:' + str(0*100/50)
    plt.text(20,7.5,txt)

def plot_animation ():
    # plot animation
    ani = animation.FuncAnimation(fig, animate, interval = 100, frames = nframe, init_func = init)
    #ani.save('anim.gif', writer="imagemagick")
    ani.save('./fig/streamline.mp4', writer="ffmpeg")
    plt.show()

def plot_image (step):
    print('Plot image:',step)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.xlim(-10,15)
    plt.ylim(-10*yscale-2.0,10*yscale+2.0)
    plt.plot([-0.5,-0.5],[-0.5,0.5],color="black",linewidth=2.5)
    plt.plot([-0.5,0.5],[0.5,0.5],color="black",linewidth=2.5)
    plt.plot([0.5,0.5],[0.5,-0.5],color="black",linewidth=2.5)
    plt.plot([0.5,-0.5],[-0.5,-0.5],color="black",linewidth=2.5)


    for line in range(0,nline):
        filename = dir_streamline + str(step) + "_" + str(line) + ".txt"
        data = np.loadtxt(filename)
        x = data[:,0]
        y = data[:,1]
        plt.plot(x,y,color = "black",linestyle="dashed")
        print("step,line:",step," ",line)

    txt = 'Steps:' + str(step) + ', ' + 'Time:' + str(0*100/50)
    plt.text(20,7.5,txt)
    filename_save = "./fig/streamline_" + str(step) + ".png"
    plt.savefig(filename_save)
    plt.show()

if __name__ == "__main__":
    #plot_animation()
    plot_image(1500)

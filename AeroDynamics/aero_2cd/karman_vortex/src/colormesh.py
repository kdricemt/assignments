%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import animation
from matplotlib.mlab import griddata

plt.rcParams["font.size"] = 30
fig = plt.figure(figsize=(40,24))
nframe = 40

def grid(x, y, z, resX=4000, resY=1000):
    "Convert 3 column data to matplotlib grid"
    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi,interp='linear')
    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z

def init():
    filename_pmap = './datafile/pmap/pmap10.txt'
    data = np.loadtxt(filename_pmap)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    X, Y, Z = grid(x, y, z)
    # 等高線図の生成。
    cont=plt.contour(X,Y,Z,48, vmin=-1.2,vmax=1.2, colors=['black'])
    #cont.clabel(fmt='%1.1f', fontsize=14)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.pcolormesh(X,Y,Z,cmap='Spectral') #カラー等高線図
    pp = plt.colorbar (orientation="vertical") # カラーバーの表示
    pp.set_label("Pressure")
    txt = 'Steps:' + str(10) + ', ' + 'Time:' + str(0*100/50)
    plt.text(20,7.5,txt)
    print("n:",10)


def animate(n):
    plt.cla() #delete previent data
    n100 = (n+1)*100
    filename_pmap = './datafile/pmap/pmap' + str(n100) + '.txt'
    data = np.loadtxt(filename_pmap)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]

    X, Y, Z = grid(x, y, z)
    # 等高線図の生成。
    cont=plt.contour(X,Y,Z,48, vmin=-1.2,vmax=1.2, colors=['black'])
    #cont.clabel(fmt='%1.1f', fontsize=14)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.pcolormesh(X,Y,Z,cmap='Spectral_r') #カラー等高線図
    txt = 'Steps:' + str(n100) + ', ' + 'Time:' + str(n100*100/50)
    plt.text(20,7.5,txt)
    print("n:",n100)

def colormap(step):
    filename_pmap = './datafile/pmap/pmap' + str(step) + '.txt'
    data = np.loadtxt(filename_pmap)
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    X, Y, Z = grid(x, y, z)
    # 等高線図の生成。
    cont=plt.contour(X,Y,Z,48, vmin=-1.2,vmax=1.2, colors=['black'])
    #cont.clabel(fmt='%1.1f', fontsize=14)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.pcolormesh(X,Y,Z,cmap='Spectral_r') #カラー等高線図
    pp = plt.colorbar (orientation="vertical") # カラーバーの表示
    pp.set_label("Pressure")
    txt = 'Time:' + str(step*100/50)
    plt.text(20,7.5,txt)
    print("n:",step)
    filename_save = "./fig/pmap" + str(step) + ".png"
    plt.savefig(filename_save)


if __name__ == "__main__":
    plotmode = "image" #image,animation
    step = 2000
    if plotmode == "image":
        colormap(step)
    else:
        ani = animation.FuncAnimation(fig, animate, interval = 100, frames = nframe, init_func = init)
        #ani.save('anim.gif', writer="imagemagick")
        ani.save('./fig/anim.mp4', writer="ffmpeg")
        plt.show()

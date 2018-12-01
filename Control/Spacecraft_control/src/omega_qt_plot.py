%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams["font.size"] = 30
fig = plt.figure(figsize=(20,15))

def plot_omega():
    plt.xlim(0,100)
    plt.ylim(-0.25,2.00)
    for n in range(1,4):
        filename_omega = './datafile/omega' + str(n) + '.txt'
        data = np.loadtxt(filename_omega)
        t = data[:,0]
        omega = data[:,1]
        #cont.clabel(fmt='%1.1f', fontsize=14)
        plt.xlabel('t')
        plt.ylabel('$\omega$[]')
        if n == 1:
            lab = "$\omega_x$"
        elif n==2:
            lab = "$\omega_y$"
        else:
            lab = "$\omega_z$"
        plt.plot(t,omega,label=lab)
        plt.legend()
    filename_save = "./fig/omega.png"
    plt.savefig(filename_save)
    plt.show()

def plot_qt():
    plt.xlim(0,100)
    plt.ylim(-1,1)
    for n in range(0,4):
        filename_qt = './datafile/qt' + str(n) + '.txt'
        data = np.loadtxt(filename_qt)
        t = data[:,0]
        qt = data[:,1]
        plt.xlabel('t')
        plt.ylabel('q[]')
        if n == 1:
            lab = "$q_0$"
        elif n==2:
            lab = "$q_1$"
        elif n==3:
            lab = "$q_2$"
        else:
            lab = "$q_3$"
        plt.plot(t,qt,label=lab)
        plt.legend()
    filename_save = "./fig/quartanion.png"
    plt.savefig(filename_save)
    plt.show()

if __name__ == "__main__":
    plot_omega()
    plot_qt()

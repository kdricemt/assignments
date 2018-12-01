%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d


def plot_omega(figsize,fsize,report_num):
    fig1 = plt.figure(figsize=figsize)
    fig2 = plt.figure(figsize=figsize)
    fig3 = plt.figure(figsize=figsize)

    if(report_num==1):
        fig = [fig1,fig1,fig1]
    else:
        fig = [fig1,fig2,fig3]

    ax1 = fig[0].add_subplot(111)
    ax2 = fig[1].add_subplot(111)
    ax3 = fig[2].add_subplot(111)

    ax = [ax1,ax2,ax3]

    for n in range(3):
        ax[n].set_xlim(0,50)
        if (n== 0):
            ax[n].set_ylim(-0.25,0.25)
        elif(n == 1):
            ax[n].set_ylim(1.85,1.925)
        elif(n == 2):
            ax[n].set_ylim(-0.25,0.25)
        if (report_num == 1):
            ax[n].set_ylim(-0.25,2)

        ax[n].set_xlabel('t')
        ax[n].set_ylabel('$\omega$[rad/s]')

        if(report_num  == 1):
            filename = './datafile/R1omega.txt'
            data = np.loadtxt(filename)
            t = data[:,0]
            omega = data[:,n+1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            if (n == 0):
                lab_1 = "$\omega_x$"
            elif(n == 1):
                lab_1 = "$\omega_y$"
            elif(n == 2):
                lab_1 = "$\omega_z$"

            ax[n].plot(t,omega,label=lab_1)
            ax[n].legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=18)
            filename_save = "./fig/R1omega.png"
            fig[n].savefig(filename_save)


        elif(report_num == 2):
            #True
            filename_true = './datafile/R2Trueomega_error.txt'
            data = np.loadtxt(filename_true)
            t = data[:,0]
            omega = data[:,n+1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            if (n== 0):
                lab_t = "$\omega_x True$"
            elif(n == 1):
                lab_t = "$\omega_y True$"
            elif(n==2):
                lab_t = "$\omega_z True$"
            ax[n].plot(t,omega,label=lab_t,color="darkviolet")

            #Est
            filename_est = './datafile/R2Estomega_error.txt'
            data = np.loadtxt(filename_est)
            t = data[:,0]
            omega = data[:,n+1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            if (n== 0):
                lab_e = "$\omega_x Est$"
            elif(n == 1):
                lab_e = "$\omega_y Est$"
            elif(n == 2):
                lab_e = "$\omega_z Est$"
            ax[n].plot(t,omega,label=lab_e,color="skyblue")

            ax[n].legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=fsize)
            filename_save = "./fig/R2omega" + str(n+1) + ".png"
            fig[n].savefig(filename_save)

def plot_qt(figsize,fsize,report_num):
    fig1 = plt.figure(figsize=figsize)
    fig2 = plt.figure(figsize=figsize)
    fig3 = plt.figure(figsize=figsize)
    fig4 = plt.figure(figsize=figsize)

    if(report_num==1):
        fig = [fig1,fig1,fig1,fig1]
    else:
        fig = [fig1,fig2,fig3,fig4]

    ax1 = fig[0].add_subplot(111)
    ax2 = fig[1].add_subplot(111)
    ax3 = fig[2].add_subplot(111)
    ax4 = fig[3].add_subplot(111)

    ax = [ax1,ax2,ax3,ax4]

    for n in range(4):
        ax[n].set_xlabel('t')
        ax[n].set_ylabel('qt[]')
        ax[n].set_xlim(0,50)
        if(n == 0):
            ax[0].set_ylim(-1,1)
        elif(n == 1):
            ax[1].set_ylim(-0.1,0.1)
        elif(n == 2):
            ax[2].set_ylim(-1,1)
        elif(n == 3):
            ax[3].set_ylim(-0.15,0.15)
        if(report_num==1):
            ax[n].set_ylim(-1,1)

        if(report_num==1):
            #Est
            filename_est = './datafile/R1qt.txt'
            data = np.loadtxt(filename_est)
            t = data[:,0]
            qt = data[:,n+1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            filename_save = "./fig/R1qt.png"

            if (n== 0):
                lab_1 = "$qt_0$"
            elif(n == 1):
                lab_1 = "$qt_1$"
            elif(n==2):
                lab_1 = "$qt_2$"
            elif(n==3):
                lab_1 = "$qt_3$"

            ax[n].plot(t,qt,label=lab_1)
            ax[n].legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=18)
            fig[n].savefig(filename_save)


        elif(report_num == 2):
            #True
            filename_true = './datafile/R2Trueqt_error.txt'
            data = np.loadtxt(filename_true)
            t = data[:,0]
            qt = data[:,n+1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            if (n== 0):
                lab_t = "$qt_0 True$"
            elif(n == 1):
                lab_t = "$qt_1 True$"
            elif(n==2):
                lab_t = "$qt_2 True$"
            elif(n==3):
                lab_t = "$qt_3 True$"

            ax[n].plot(t,qt,label=lab_t,color="darkviolet")

            #Est
            filename_est = './datafile/R2Estqt_error.txt'
            data = np.loadtxt(filename_est)
            t = data[:,0]
            qt = data[:,n+1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            filename_save = "./fig/R2qt" + str(n) + ".png"

            if (n== 0):
                lab_e = "$qt_0 Est$"
            elif(n == 1):
                lab_e = "$qt_1 Est$"
            elif(n==2):
                lab_e = "$qt_2 Est$"
            elif(n==3):
                lab_e = "$qt_3 Est$"

            ax[n].plot(t,qt,label=lab_e,color="skyblue")
            ax[n].legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=fsize)
            fig[n].savefig(filename_save)

def plot_p(figure,fsize):
        fig1 = plt.figure(figsize=figsize)
        fig2= plt.figure(figsize=figsize)
        fig3 = plt.figure(figsize=figsize)
        fig4 = plt.figure(figsize=figsize)
        fig5 = plt.figure(figsize=figsize)
        fig6 = plt.figure(figsize=figsize)
        fig7= plt.figure(figsize=figsize)

        ax1 = fig1.add_subplot(111)
        ax2 = fig2.add_subplot(111)
        ax3 = fig3.add_subplot(111)
        ax4= fig4.add_subplot(111)
        ax5 = fig5.add_subplot(111)
        ax6 = fig6.add_subplot(111)
        ax7 = fig7.add_subplot(111)


        ax = [ax1,ax2,ax3,ax4,ax5,ax6,ax7]
        fig = [fig1,fig2,fig3,fig4,fig5,fig6,fig7]

        for n in range(7):
            ax[n].set_xlabel('t')
            ax[n].set_ylabel('P - Estimate error covariance')
            ax[n].set_xlim(0,50)
            ax[n].set_ylim(-0.1,0.1)

            filename_p = './datafile/R2p.txt'
            data_p = np.loadtxt(filename_p)
            t = data_p[:,0]
            p = data_p[:,n+1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            lab_p = "$\pm \sqrt{p}$" + str(n) + str(n) ;

            filename_e = './datafile/R2e.txt'
            data_e = np.loadtxt(filename_e)
            e = data_e[:,n+1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            if n < 4:
                lab_e = "$\delta$p" + str(n);
            elif n == 4:
                lab_e = "$\delta \omega_x$"
            elif n == 5:
                lab_e = "$\delta \omega_y$"
            elif n == 6:
                lab_e = "$\delta \omega_z$"

            ax[n].plot(t,p,label=lab_p,color="darkviolet")
            ax[n].plot(t,e,label=lab_e,color="skyblue")
            ax[n].plot(t,-p,color="darkviolet")

            filename_save = "./fig/R2pmatrix" + str(n) + ".png"

            ax[n].legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=fsize)
            fig[n].savefig(filename_save)

def plot_k(figure,fsize):
        fig1 = plt.figure(figsize=figsize)

        ax1 = fig1.add_subplot(111)
        ax2 = fig1.add_subplot(111)
        ax3 = fig1.add_subplot(111)

        ax = [ax1,ax2,ax3]

        for n in range(3):

            ax[n].set_xlabel('t')
            ax[n].set_ylabel('K - Kalman Gain')
            ax[n].set_xlim(0,100)
            ax[n].set_ylim(0,2.0)

            #True
            filename_k = './datafile/R2k_' + str(n) + '.txt'
            data = np.loadtxt(filename_k)
            t = data[:,0]
            k = data[:,1]
            #cont.clabel(fmt='%1.1f', fontsize=14)
            if n==0:
                lab_k = "$\sigma_v = \sigma_w = 0.02$";
            elif n==1:
                lab_k = "$\sigma_v(0.03) > \sigma_w(0.01)$";
            else:
                lab_k = "$\sigma_v(0.01) < \sigma_w(0.03)$";

            ax[n].plot(t,k,label=lab_k)

            filename_save = "./fig/R2kmatrix_compare.png"

            ax[n].legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=fsize)

            fig1.savefig(filename_save)

def plot_dcm(figsize,fsize):
    fig1 = plt.figure(figsize=figsize)
    fig2= plt.figure(figsize=figsize)
    fig3 = plt.figure(figsize=figsize)
    ax1 = fig1.gca(projection='3d')
    ax2 = fig2.gca(projection='3d')
    ax3 = fig3.gca(projection='3d')
    fig = [fig1,fig2,fig3]
    ax = [ax1,ax2,ax3]

    data1 = np.loadtxt("./datafile/R1dcm_column0.txt")
    data2 = np.loadtxt("./datafile/R1dcm_column1.txt")
    data3 = np.loadtxt("./datafile/R1dcm_column2.txt")

    t = data1[:,0]
    col0_0 = data1[:,1]
    col0_1 = data1[:,2]
    col0_2 = data1[:,3]
    col1_0 = data2[:,1]
    col1_1 = data2[:,2]
    col1_2 = data2[:,3]
    col2_0 = data3[:,1]
    col2_1 = data3[:,2]
    col2_2 = data3[:,3]

    base_1 = np.array([1,0,0])
    base_2 = np.array([0,1,0])
    base_3 = np.array([0,0,1])
    v1_0 = []
    v1_1 = []
    v1_2 = []
    v2_0= []
    v2_1= []
    v2_2= []
    v3_0 = []
    v3_1= []
    v3_2= []

    filename = './datafile/R1omega.txt'
    data = np.loadtxt(filename)
    t = data[:,0]
    omega_x = data[:,1]
    omega_y = data[:,2]
    omega_z = data[:,3]
    I = np.array([1.9,1.6,2.0])
    omega_0 = np.array([omega_x[0],omega_y[0],omega_z[0]])
    H = I * omega_0 #angular momentum
    H = H * 5/ np.linalg.norm(H)
    print(H)

    for i in range(len(t)):
        if((i%5) == 0):
            dcm = np.array([[col0_0[i],col1_0[i],col2_0[i]],
                          [col0_1[i],col1_1[i],col2_1[i]],
                          [col0_2[i],col1_2[i],col2_2[i]]])
            v1 = np.dot(dcm,base_1);
            v2 = np.dot(dcm,base_2);
            v3 = np.dot(dcm,base_3);
            v1_0.append(v1[0])
            v1_1.append(v1[1])
            v1_2.append(v1[2])

            v2_0.append(v2[0])
            v2_1.append(v2[1])
            v2_2.append(v2[2])

            v3_0.append(v3[0])
            v3_1.append(v3[1])
            v3_2.append(v3[2])

    ax1.set_xlabel("x")
    ax1.set_ylabel("y")
    ax1.set_zlabel("z")
    ax2.set_xlabel("x")
    ax2.set_ylabel("y")
    ax2.set_zlabel("z")
    ax3.set_xlabel("x")
    ax3.set_ylabel("y")
    ax3.set_zlabel("z")

    ax1.plot(v1_0,v1_1,v1_2,color='darkorchid',label="x_axis")
    ax2.plot(v2_0,v2_1,v2_2,color='darkorchid',label="y_axis")
    ax3.plot(v3_0,v3_1,v3_2,color='darkorchid',label="z_axis")

    ax1.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=fsize)
    ax2.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=fsize)
    ax3.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=fsize)

    ax1.view_init(30,220)
    ax2.view_init(30,220)
    ax3.view_init(30,220)

    fig1.savefig("./fig/R1_xaxis.png")
    fig2.savefig("./fig/R1_yaxis.png")
    fig3.savefig("./fig/R1_zaxis.png")


if __name__ == "__main__":
    fsize = 12 #font size of legend
    plt.rcParams["font.size"] = 10
    figsize = [12,6]
    # for report 1
    #plot_omega(figsize,fsize,1)
    #plot_qt(figsize,fsize,1)
    #plot_dcm([8,8],fsize)
    # for report 2
    #plot_omega(figsize,fsize,2)
    #plot_qt(figsize,fsize,2)
    #plot_p(figsize,fsize)
    plot_k(figsize,fsize)

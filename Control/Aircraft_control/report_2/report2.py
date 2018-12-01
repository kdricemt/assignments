%matplotlib inline
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import sympy as sym

#parameter
m = 7500 #[kg]
s = 27.9 #[m^2]
cl_alpha = 4.30
cd_0 = 0.0548
k = 3.02
rho = 0.74 #[kg/m^3]
g = 9.806 #[m_s^2]

alpha = 0 #radian
thrust = 0

#結果格納用array
t_array = []
u_array = []
gamma_array = []
psi_array = []
phi_array = []
xe_array = []
ye_array = []
ze_array = []

pi = 3.141592
rad = pi/180

def cal_alpha_thrust():
    global m,s,rho,g,cl_alpha,cd_0,k,alpha,thrust,pi,rad
    rad = pi / 180

    (alpha_s,thrust_s) = sym.symbols("alpha_s thrust_s")
    u = 180
    q = 0.5 * rho * u**2

    f1 = - q*s*(cd_0 + k*alpha_s**2) + thrust_s * (1 - alpha_s**2/2)
    f2 = q * s * cl_alpha * alpha_s + thrust_s*alpha_s - m*g
    res = sym.solve([f1,f2],[alpha_s,thrust_s])
    alpha = res[0][0]
    thrust = res[0][1]
    print("alpha=",alpha/rad,"[deg]")
    print("Thrust=",thrust,"[N]")

#ロール角
def phi(time):
    return 0

#速度
def diff_u(t,u,gamma,psi):
    global m,s,rho,g,cl_alpha,cd_0,k,alpha,thrust,pi,rad
    q = 0.5 * rho * u**2
    drag = q * s * (cd_0 + k*alpha**2)
    diff_u = (-drag + thrust*sym.cos(alpha))/m - g * sym.sin(gamma)
    return diff_u

#経路角
def diff_gamma(t,u,gamma,psi):
    global m,s,rho,g,cl_alpha,cd_0,k,alpha,thrust,pi,rad
    q = 0.5 * rho * u**2
    lift = q * s * cl_alpha * alpha
    diff_gamma = ((lift + thrust*sym.sin(alpha))* sym.cos(phi(t))/m - g*sym.cos(gamma)) / u
    return diff_gamma

#yaw
def diff_psi(t,u,gamma,psi):
    global m,s,rho,g,cl_alpha,cd_0,k,alpha,thrust,pi,rad
    q = 0.5 * rho * u**2
    lift = q * s * cl_alpha * alpha
    diff_psi = ((lift + thrust*sym.sin(alpha))/m) * (sym.sin(phi(t))/sym.cos(gamma)) / u
    return diff_psi

def diff_xe(t,u,gamma,psi):
    diff_xe = u * sym.cos(gamma) * sym.cos(psi)
    return diff_xe

def diff_ye(t,u,gamma,psi):
    diff_ye = u * sym.cos(gamma) * sym.sin(psi)
    return diff_ye

def diff_ze(t,u,gamma,psi):
    diff_ze = -u * sym.sin(gamma)
    return diff_ze

def cal_flightroute():
    global m,s,u,rho,g,cl_alpha,cd_0,k,alpha,thrust,pi,rad
    global xe_array,ye_array,ze_array,t_array,u_array,gamma_array,psi_array,phi_array

    #初期化
    gamma = 0 #経路角
    psi = 0 #ヨー角
    phi_data=0 #ロール
    t=0
    dt = 0.5 #刻み幅
    u=200
    endt = 400.0 #終了時刻
    xe = 0
    ye = 0
    ze = -5000
    index=1

    #結果を格納する配列
    t_array = [t]
    u_array = [u]
    gamma_array = [gamma]
    psi_array = [psi]
    phi_array = [phi_data]
    xe_array = [xe]
    ye_array = [ye]
    ze_array = [-ze]
    #周期測定用
    index_array_gamma = []
    index_array_u = []

    #簡単のため名前を置き換える
    x = u
    y = gamma
    z = psi

    #ルンゲクッタ計算用array
    k0 = [0,0,0,0,0,0]
    k1 = [0,0,0,0,0,0]
    k2 = [0,0,0,0,0,0]
    k3 = [0,0,0,0,0,0]

    while t < endt:
        k0[0]=dt*diff_u(t,x,y,z)
        k0[1]=dt*diff_gamma(t,x,y,z)
        k0[2]=dt*diff_psi(t,x,y,z)
        k0[3]=dt*diff_xe(t,x,y,z)
        k0[4]=dt*diff_ye(t,x,y,z)
        k0[5]=dt*diff_ze(t,x,y,z)

        k1[0]=dt*diff_u(t+dt/2.0,x+k0[0]/2.0,y+k0[1]/2.0,z+k0[2]/2.0)
        k1[1]=dt*diff_gamma(t+dt/2.0,x+k0[0]/2.0,y+k0[1]/2.0,z+k0[2]/2.0)
        k1[2]=dt*diff_psi(t+dt/2.0,x+k0[0]/2.0,y+k0[1]/2.0,z+k0[2]/2.0)
        k1[3]=dt*diff_xe(t+dt/2.0,x+k0[0]/2.0,y+k0[1]/2.0,z+k0[2]/2.0)
        k1[4]=dt*diff_ye(t+dt/2.0,x+k0[0]/2.0,y+k0[1]/2.0,z+k0[2]/2.0)
        k1[5]=dt*diff_ze(t+dt/2.0,x+k0[0]/2.0,y+k0[1]/2.0,z+k0[2]/2.0)

        k2[0]=dt*diff_u(t+dt/2.0,x+k1[0]/2.0,y+k1[1]/2.0,z+k1[2]/2.0)
        k2[1]=dt*diff_gamma(t+dt/2.0,x+k1[0]/2.0,y+k1[1]/2.0,z+k1[2]/2.0)
        k2[2]=dt*diff_psi(t+dt/2.0,x+k1[0]/2.0,y+k1[1]/2.0,z+k1[2]/2.0)
        k2[3]=dt*diff_xe(t+dt/2.0,x+k1[0]/2.0,y+k1[1]/2.0,z+k1[2]/2.0)
        k2[4]=dt*diff_ye(t+dt/2.0,x+k1[0]/2.0,y+k1[1]/2.0,z+k1[2]/2.0)
        k2[5]=dt*diff_ze(t+dt/2.0,x+k1[0]/2.0,y+k1[1]/2.0,z+k1[2]/2.0)

        k3[0]=dt*diff_u(t+dt,x+k2[0],y+k2[1],z+k2[2])
        k3[1]=dt*diff_gamma(t+dt,x+k2[0],y+k2[1],z+k2[2])
        k3[2]=dt*diff_psi(t+dt,x+k2[0],y+k2[1],z+k2[2])
        k3[3]=dt*diff_xe(t+dt,x+k2[0],y+k2[1],z+k2[2])
        k3[4]=dt*diff_ye(t+dt,x+k2[0],y+k2[1],z+k2[2])
        k3[5]=dt*diff_ze(t+dt,x+k2[0],y+k2[1],z+k2[2])

        x=x+(k0[0]+2.0*k1[0]+2.0*k2[0]+k3[0])/6.0
        y=y+(k0[1]+2.0*k1[1]+2.0*k2[1]+k3[1])/6.0
        z=z+(k0[2]+2.0*k1[2]+2.0*k2[2]+k3[2])/6.0
        xe=xe+(k0[3]+2.0*k1[3]+2.0*k2[3]+k3[3])/6.0
        ye=ye+(k0[4]+2.0*k1[4]+2.0*k2[4]+k3[4])/6.0
        ze=ze+(k0[5]+2.0*k1[5]+2.0*k2[5]+k3[5])/6.0

        t_array.append(t)
        u_array.append(x)
        gamma_array.append(y/rad)
        psi_array.append(z/rad)
        phi_data = phi(t)
        phi_array.append(phi_data/rad)
        xe_array.append(xe)
        ye_array.append(ye)
        ze_array.append(-ze)

        if index >0:
            if u_array[index]<u_array[index-1] and u_array[index-1]>u_array[index-2]:
                index_array_u.append(index)
            if gamma_array[index]<gamma_array[index-1] and gamma_array[index-1]>gamma_array[index-2]:
                index_array_gamma.append(index)

        t += dt
        index += 1

    print("U=",u_array)
    print("x_e:",xe_array)
    print("y_e:",ye_array)
    print("z_e:",ze_array)
    print("gamma:",gamma_array)
    print("psi:",psi_array)
    print("phi:",phi_array)

    print("t_gamma:",index_array_gamma)
    print("t_u:",index_array_u)
    print("gamma_max:")
    for index_g in index_array_gamma:
        print(gamma_array[index_g])
    print("u_max:")
    for index_u in index_array_u:
        print(u_array[index_u])

def plot_route():
    global t_array,u_array,xe_array,ye_array,ze_array,gamma_array,psi_array,phi_array
    fig2 = plt.figure()
    fig6 = plt.figure()

    plt.rcParams["font.size"] = 10
    #create object
    ax2 = fig2.add_subplot(111)
    ax6 = fig6.add_subplot(111)

    # Data for a 3D line
    ax2.plot(t_array,u_array,color='red')
    ax6.plot(t_array,gamma_array,color='red')

    # set title
    ax2.set_title("U")
    ax6.set_title("gamma[deg]")

    #set label
    ax2.set_xlabel("t[s]")
    ax2.set_ylabel("U[m/s]")
    ax6.set_xlabel("t[s]")
    ax6.set_ylabel("gamma[deg]")

    filename2 = "U_2.png"
    fig2.savefig(filename2)
    filename6 = "gamma_2.png"
    fig6.savefig(filename6)

if __name__ == '__main__':
    print("Q1")
    cal_alpha_thrust()
    print("Q2")
    cal_flightroute()
    plot_route()

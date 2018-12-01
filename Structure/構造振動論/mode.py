%matplotlib inline
import numpy as np
import math
import matplotlib.pyplot as plt
import sympy as sym
import random
import scipy.linalg

rho = 2.7 * 10**(3) #[kg/m^3]
e = 70 * 10**(9) #Pa
h = 0.001 #[m]
l = 0.1 #[m]
b = 1 #[m]
name = "simpson"

gamma = [4.73004,7.85320,10.9956,14.1371]
alpha = [0,0,0,0]


def cal_alpha():
    global alpha
    for i in range(0,4):
        g = gamma[i]
        alpha[i] = (math.cosh(g) - math.cos(g))/(math.sinh(g)-math.sin(g))
    print("alpha")
    print(alpha)

#n is mode number(0,1,2)
def phi(n,xi):
    n = int(n)
    global gamma,alpha
    gamma_gsi = gamma[n]*xi
    ans = math.cosh(gamma_gsi) + math.cos(gamma_gsi) - alpha[n]*(math.sinh(gamma_gsi) + math.sin(gamma_gsi))
    return ans

def a_function(xi,i,j):
    global rho,e,h,l,b
    return rho * b * h * (1 + 2 * xi) * phi(i,xi) * phi(j,xi)

def a_integral(i,j):
    a_integ = simpson(a_function,i,j,0,1)
    return a_integ

def phi_dash(n,xi):
    n = int(n)
    global gamma,alpha
    gamma_gsi = gamma[n]*xi
    ans = math.cosh(gamma_gsi) - math.cos(gamma_gsi) - alpha[n]*(math.sinh(gamma_gsi) - math.sin(gamma_gsi))
    return ans

def b_function(xi,i,j):
    global rho,e,h,l,b
    return  gamma[i]**2 * gamma[j]**2 * e * b * h**3 / (12 * l**4) * (1 + 2 * xi)**3 * phi_dash(i,xi) * phi_dash(j,xi)

def b_integral(i,j):
    b_integ = simpson(b_function,i,j,0,1)
    return b_integ

def omega(n,x):
    global l,gamma,alpha
    gamma_x = gamma[n] * x / l
    omega = (math.cosh(gamma_x) + math.cos(gamma_x) - alpha[n]*(math.sinh(gamma_x) + math.sin(gamma_x)))/math.sqrt(l)
    return omega

def mode(x,vr):
    sum_mode = 0
    for i in range(0,4):
        sum_mode = sum_mode + vr[0][i] * omega(i,x)
    return sum_mode

def lamda_uni(n):
    global rho,e,h,l,b,gamma
    i = b * (2*h)**3 /12
    a = b * (2*h)
    return (gamma[n]/l)**2 * math.sqrt(e*i/(rho*a))

def solve():
    global l
    A = np.zeros((4,4))
    B = np.zeros((4,4))
    for i in range(0,4):
        for j in range(0,4):
            A[i,j] = a_integral(i,j)
            B[i,j] = b_integral(i,j)
    print('')
    print('A')
    print("[",A[0,0]," ",A[0,1]," ",A[0,2],"]")
    print("[",A[1,0]," ",A[1,1]," ",A[1,2],"]")
    print("[",A[2,0]," ",A[2,1]," ",A[2,2],"]")
    print('B')
    print("[",B[0,0]," ",B[0,1]," ",B[0,2],"]")
    print("[",B[1,0]," ",B[1,1]," ",B[1,2],"]")
    print("[",B[2,0]," ",B[2,1]," ",B[2,2],"]")
    print("")

    w, vr = scipy.linalg.eig(B, A)
    lamda = np.sort(np.sqrt(w))
    #f = omega/(2*np.pi)
    for i in range(len(vr)): #正規化
        vr [i] = vr[i]/np.linalg.norm(vr[i])
    print("")
    print("■ 変断面はり")
    print("固有振動数 λ:",lamda)
    print("固有ベクトル a:")
    vr1,vr2,vr3,vr4 = np.zeros(4)
    vr1 = vr[0:1]
    vr2 = vr[1:2]
    vr3 = vr[2:3]
    vr4 = vr[3:4]
    #vr3[0][0] = -vr3[0][0]
    #vr3[0][1]= -vr3[0][1]
    print("1次モード",vr1)
    print("2次モード",vr2)
    print("3次モード",vr3)
    print("4次モード",vr4)

    # λ for uniform_spansize
    print("")
    print("■ 一様断面はり")
    print("固有振動数 λ")
    for i in range(0,3):
        print("モード:",i+1,"  λ:",lamda_uni(i))

    #plot
    x_array = np.linspace(0,l,50)
    v_mode = np.vectorize(mode)

    #fig1 arbitrary spansize
    fig1 = plt.figure()
    y_array1 = np.zeros(50)
    y_array2 = np.zeros(50)
    y_array3 = np.zeros(50)
    index = 0

    for x in x_array:
        y_array1[index] = mode(x,vr4)
        y_array2[index] = mode(x,vr1)
        y_array3[index] = mode(x,vr2)
        index = index + 1

    print('')
    plt.title("arbitary spansize")
    plt.xlabel("x")
    plt.plot(x_array,y_array1,color="blue",label="Mode 1")
    plt.plot(x_array,y_array2,color="orange",label="Mode 2")
    plt.plot(x_array,y_array3,color="green",label="Mode 3")
    plt.legend()
    figname1 = 'arbitary_spansize.png'
    plt.savefig(figname1)

    #fig2 uniform_spansize
    fig2 = plt.figure()
    y_array1_uni = np.zeros(50)
    y_array2_uni = np.zeros(50)
    y_array3_uni = np.zeros(50)
    index2 = 0

    for x in x_array:
        y_array1_uni[index2] = omega(0,x)
        y_array2_uni[index2] = omega(1,x)
        y_array3_uni[index2] = omega(2,x)
        index2 = index2 + 1

    plt.title("uniform spansize")
    plt.xlabel("x")
    plt.plot(x_array,y_array1_uni,color="blue",label="Mode 1")
    plt.plot(x_array,y_array2_uni,color="orange",label="Mode 2")
    plt.plot(x_array,y_array3_uni,color="green",label="Mode 3")
    plt.legend()
    figname2 = 'uniform_spansize.png'
    plt.savefig(figname2)



#simpson's formula
def simpson(f,i,j,start,finish):
    global name
    name = "simpson"
    N = 10
    h = (finish-start)/N
    s = (h/3) * sum((f(h*k,i,j) + 4*f(h*(k+1),i,j) + f(h*(k+2),i,j)) for k in range(0,N-1, 2))
    if f == a_function:
        name = "a"
    else:
        name = "b"
    #print("simpson function:",name,"  i,j:",i,",",j,"  sum:", s)
    return s

if __name__ == "__main__":
    cal_alpha()
    solve()

#for i in range(0,5):
#    print(a_integral(2,i))
#p_matrix()

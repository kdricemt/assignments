%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt


# 1D hyperbolic型の方程式の解
def plot_1D(scheme_str):
    filepath = "./datafile/1Dhyperbolic/"
    figpath = "./fig/1Dhyperbolic/"
    filename = filepath + scheme_str + ".txt"
    figname =  figpath + scheme_str + ".png"
    file_exact = "./datafile/1Dhyperbolic/exact.txt"
    # import data
    data = np.loadtxt(filename)
    x = data[:,0]
    y = data[:,1]
    data_exact = np.loadtxt(file_exact)
    y_exact = data_exact[:,1]
    # plot data
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x,y,color="black", marker="o", markersize=5, markeredgewidth=3, markeredgecolor="black",
  markerfacecolor="black",label="Numerical Solution")
    ax.plot(x,y_exact, color="black",linestyle="dashed", marker="o", markersize=5, markeredgewidth=3, markeredgecolor="black",
  markerfacecolor="black",label="Exact Solution")
    ax.set_title(scheme_str)
    ax.set_xlim([0.4,1])
    ax.legend()
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.savefig(figname)
    return fig,ax

def plot_tvd(scheme_str):
    filepath = "./datafile/sym_tvd/"
    figpath = "./fig/sym_tvd/"
    filename = filepath + scheme_str + ".txt"
    figname =  figpath + scheme_str + ".png"
    file_exact = "./datafile/sym_tvd/exact.txt"
    # import data
    data = np.loadtxt(filename)
    x = data[:,0]
    y = data[:,1]
    data_exact = np.loadtxt(file_exact)
    y_exact = data_exact[:,1]
    # plot data
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.plot(x,y,color="black", marker="o", markersize=5, markeredgewidth=3, markeredgecolor="black",
  markerfacecolor="black",label="Numerical Solution")
    ax.plot(x,y_exact, color="black",linestyle="dashed", marker="o", markersize=5, markeredgewidth=3, markeredgecolor="black",
  markerfacecolor="black",label="Exact Solution")
    ax.set_title(scheme_str)
    ax.legend()
    ax.set_xlim([0.4,1])
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.savefig(figname)
    return fig,ax


scheme_str_array = ["center","upper_stream","Lax_Wendroff"]
for scheme_1d in scheme_str_array:
    plot_1D(scheme_1d)
tvd_scheme_str_array = ["minmod1_ecp0.1","minmod1_ecp0.2","minmod1_ecp0.3","minmod2_ecp0.1","superbee_ecp0.1"]
for scheme_tvd in tvd_scheme_str_array:
    plot_tvd(scheme_tvd)

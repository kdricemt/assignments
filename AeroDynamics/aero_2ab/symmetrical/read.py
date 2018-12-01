#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.mlab import griddata

matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
plt.rcParams["font.size"] = 35
#plt.axes().set_aspect('equal')
fig = plt.figure(figsize=(50,30))

def grid(x, y, z, resX=1000, resY=1000):
    "Convert 3 column data to matplotlib grid"
    xi = np.linspace(min(x), max(x), resX)
    yi = np.linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi,interp='linear')
    X, Y = np.meshgrid(xi, yi)
    return X, Y, Z

def wing_shape1():
    data = np.loadtxt('xy_zeta.txt',delimiter='\t')
    x = data[:,0]
    y = data[:,1]
    plt.plot(x,y,color='black',linestyle='solid')

def wing_shape2():
    data = np.loadtxt('xy_zeta.txt',delimiter='\t')
    x = data[:,0]
    y = data[:,1]
    plt.plot(x,y,color='black',linestyle='solid')

def streamline():
    data = np.loadtxt('streamline.txt',delimiter='\t')
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    X, Y, Z = grid(x, y, z)
    levels = np.arange(-2.0,2.0,0.1)
    CF = plt.contour(X, Y, Z,levels,linetypes=8,cmap = cm.spring)
    CB = plt.colorbar(CF)
    CB.set_label('f_image') 
    plt.xlim(-2, 2)
    plt.xlabel("x")
    plt.ylim(-1.5, 1.5)
    
def pressure():
    data = np.loadtxt('cp.txt',delimiter='\t')
    x = data[:,0]
    y = data[:,1]
    z = data[:,2]
    X, Y, Z = grid(x, y, z)
    levels = np.arange(-2.0,2.0,0.1)
    CF = plt.contour(X, Y, Z,levels,linetypes=8,cmap = cm.autumn)
    CB = plt.colorbar(CF)
    CB.set_label('Cp') 
    plt.xlim(-2, 2)
    plt.xlabel("x")
    plt.ylim(-1.5, 1.5)

def display(graph):   
    if(graph==0):
        wing_shape1()
        streamline()
        filename = "stream_j.png"
        plt.savefig(filename)
    else:
        wing_shape2()
        pressure()
        filename = "pressure_j.png"
        plt.savefig(filename)
         
display(1)

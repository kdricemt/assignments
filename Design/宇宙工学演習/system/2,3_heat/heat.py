#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import matplotlib.pyplot as plt

p_s = 1358
alpha_s = 0.2
rad = math.pi/180

#stephan_boltzman_constant
sigma = 5.67 * math.pow(10,-8)
temp_wall = 20 + 273.15
epsilon = 0.8
default_ratio = 1
pad_ratio = 0.9
tar_p_ratio = 0.9
beta_sum = 23.4 
beta_spr = 0
dalpha = 0.1 #積分用
alpha = 0

q_tarp = 0
q_tarm = 0
q_sunp = 0
q_sunm = 0
q_padp = 0
q_padm = 0

def init_parameter():
    global alpha,q_tarp,q_tarm,q_sunp,q_sunm,q_padp,q_padm
    alpha = 0
    q_tarp = 0
    q_tarm = 0
    q_sunp = 0
    q_sunm = 0
    q_padp = 0
    q_padm = 0

def cal_radiation(beta):
    global alpha,q_tarp,q_tarm,q_sunp,q_sunm,\
    q_padp,q_padm,dalpha
    
    while alpha < 360:
        q_tarp += max(0,-math.cos(alpha*rad) * dalpha)
        q_tarm += max(0,math.cos(alpha*rad) * dalpha)
        q_sunp += max(0,math.sin(alpha*rad) * dalpha)
        q_sunm += max(0, -math.sin(alpha*rad) * dalpha)
        q_padp += max(0,1 * dalpha)
        q_padm += max(0,0 * dalpha)
        alpha += dalpha
    
    #係数をかける
    q_tarp *= alpha_s * p_s * math.cos(beta*rad)
    q_tarm *= alpha_s * p_s * math.cos(beta*rad)
    q_sunp *= alpha_s * p_s * math.cos(beta*rad)
    q_sunm *= alpha_s * p_s * math.cos(beta*rad)
    q_padp *= alpha_s * p_s * math.sin(beta*rad)
    q_padm *= alpha_s * p_s * math.sin(beta*rad)
    
    #一周の平均
    q_tarp *= 1/360.0
    q_tarm *= 1/360.0
    q_sunp *= 1/360.0
    q_sunm *= 1/360.0
    q_padp *= 1/360.0
    q_padm *= 1/360.0
    
    #放射能力
    global epsilon,tar_p_ratio,pad_ratio,default_ratio
    
    c1 =  epsilon * sigma * math.pow(temp_wall,4)
    p_tarp = c1 * tar_p_ratio - q_tarp
    p_tarm = c1 * default_ratio - q_tarm
    p_sunp = c1 * default_ratio - q_sunp
    p_sunm = c1 * default_ratio - q_sunm
    p_padp = c1 * pad_ratio - q_padp
    p_padm = c1 * pad_ratio - q_padm
    
    print("beta",beta)
    
    print("q_tarp",q_tarp)
    print("q_tarm",q_tarm)
    print("q_sunp",q_sunp)
    print("q_sunm",q_sunm)
    print("q_padp",q_padp)
    print("q_padm",q_padm)
    
    print("p_tarp",p_tarp)
    print("p_tarm",p_tarm)
    print("p_sunp",p_sunp)
    print("p_sunm",p_sunm)
    print("p_padp",p_padp)
    print("p_padm",p_padm)
    print("\n")
    
if __name__ == "__main__":
    init_parameter()
    cal_radiation(beta_sum)
    init_parameter()
    cal_radiation(beta_spr)
    
    



    
        
        


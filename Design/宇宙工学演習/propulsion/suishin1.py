#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math
import sympy

density_water = 1000 #ρ　kg/m^3
diameter = 20 *math.pow(10,-3) #D
bulk_velosity = 0.05 #m/s ub
kinematic_viscocity = 0.668*math.pow(10,-6) #v 動粘性係数
prandtl_num = 4.39 #プラントル数 Pr
specific_heat = 4.217 #水の比熱c
wall_temp = 100 #壁の温度
enter_temp = 20
exit_temp = 60

def viscocity(): #粘性係数　μ
    return density_water * kinematic_viscocity

def reynolds_num(): #レイノルズ数　Re
    return bulk_velosity * diameter / kinematic_viscocity

def thermal_conductivity(): #熱伝導率 k
    k = viscocity() * specific_heat / prandtl_num
    return k

def lmtd(): #対数平均温度差　ΔTm
    delta_t1 = wall_temp - enter_temp
    delta_t2 = wall_temp - exit_temp
    lmtd = 1.0*(delta_t1-delta_t2)/math.log(1.0*delta_t1/delta_t2)
    return lmtd

def quantity_heat(): #熱量Q
    q = math.pi * math.pow(diameter,2) * density_water *\
                bulk_velosity * specific_heat * (exit_temp - enter_temp) / 4
    return q

def nusselt_number(length):
    gz = diameter * reynolds_num() * prandtl_num / length       
    return 3.65 + (0.0668 * gz /(1 + 0.04 * math.pow(gz,2/3)))
    
def heat_transfer_coefficient(length): #熱伝達率
    nusselt_num = nusselt_number(length)
    h = nusselt_num * thermal_conductivity() / diameter
    return h

def solve():
    INF = 10000
    min_eq = INF
    pipe_length = INF
    length = sympy.symbols('length')
    for length in range(1,200,1):
        eq = quantity_heat() /(math.pi * diameter * (1.0*length/10) * lmtd())\
                -heat_transfer_coefficient(1.0*length/10)
        if abs(eq)<min_eq:
            min_eq = abs(eq)
            pipe_length = length*1.0/10
        else:
            min_eq = min_eq
    print("tube length:",pipe_length)
    print("diff;",min_eq)
    print("nusselt:",nusselt_number(pipe_length))
    print("imtd:",lmtd())
    print("Reynolds_num:",reynolds_num())
    
if __name__ == '__main__':
    solve()


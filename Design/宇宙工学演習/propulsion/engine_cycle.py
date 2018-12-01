#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math
import decimal
mpl.rcParams['font.family'] = 'AppleGothic'
# 設定,要求したパラメータ

chamber_pressure = 20.64 * math.pow(10, 6)  # N/m^2,SSME 燃焼室圧力 cp
Pa = 1013 * math.pow(10, 2)  # 大気圧
P0 = 20.64 * math.pow(10, 6)  # SSME Pressure inlet
PaP0ratio = Pa / P0  # 最適膨張を仮定
delta_v = 7000  # ΔV
ps = 0.3  # 標準混合推進剤密度
finerts = 0.1  # 標準構造質量比
payload_mass = 2000  # kg
ts = 0.01  # 燃焼室滞在時間 10ms
a_i = 1.3 * 9.814  # 初期加速度 1.3G
eta_j = 0.96  # ノズル効率

# 定数
gamma = 1.4  # 比熱比γ
g = 9.8  # 重力加速度
R = 8.314 * 1000  # 気体定数
H2_density = 0.071  # g/cm3 liquid
O2_density = 1.14  # g/cm^3 liquid
F2_density = 1.51  # g/cm^3 liquid
cp = 8.31 * (gamma) / ((gamma) - 1)  # 定圧比熱

#mass_ratio = MR
# adiabatic_Temp = Tf 断熱火炎温度
# Vj 排気速度
# density = ρ 混合推進剤密度
# finert 構造質量比
# thrust = F 推力
# isp 比推力
# thrust_coefficient = Cf 推力係数
# a_i 初期加速度 1.3G
# throatA = A*,At スロート断面積

def mr_Change(mixture_ratio):
    # 質量混合比をモルに変える moxider/mfuel
    moll_ratio = mixture_ratio / 16
        # Oの原子量
    return moll_ratio

# 全てH21mol基準で


def product_H2_0(moll_ratio):  # H2がどれぐらい残るのか
    # H2+0.5O2=H20
    if moll_ratio > 0.5:
        return 0
    else:
        return 1 - (2 * moll_ratio)  # product(生成物)のmoll数を返す


def product_O2(moll_ratio):
    if moll_ratio > 0.5:
        return moll_ratio - 0.5
    else:
        return 0


def product_H2o(moll_ratio):
    if moll_ratio > 0.5:
        return 1
    else:
        return 2 * moll_ratio

def Adiabatic_Temp(mixture_ratio,entrance_temp): #entrance_tempは主燃焼室燃焼前温度
    moll_ratio = mr_Change(mixture_ratio)
    pH2_0 = product_H2_0(moll_ratio)
    pO2 = product_O2(moll_ratio)
    pH2o = product_H2o(moll_ratio)
    enth = pH2_0 * cp + pO2 * cp + pH2o * cp
    hof = 241.9 * pH2o * 1000  # 単位J 生成熱
    temp = hof / enth + entrance_temp
    if temp > 2500:
            temp = (2 * hof + entrance_temp * enth) / (enth + 0.0004 * hof)
    return temp
   


def Products_mass(mixture_ratio):
    # 水素1molとした時の反応後の重さ/mol
    moll_ratio = mr_Change(mixture_ratio)
    pH2_0 = product_H2_0(moll_ratio)
    pO2 = product_O2(moll_ratio)
    pH2o = product_H2o(moll_ratio)
    m = pH2_0 * 2 + pO2 * 32 + pH2o * 18 / (pH2_0 + pO2 + pH2o)
   
    return m


def Vj(mixture_ratio,entrance_temp):
    m = Products_mass(mixture_ratio)
    temp = Adiabatic_Temp(mixture_ratio,entrance_temp)
    gammaratio = (gamma - 1) / gamma
    powg = math.pow(PaP0ratio, gammaratio)
    vj = math.sqrt(2 * eta_j * (gamma / (gamma - 1))
                   * (R / m) * temp * (1 - powg))
    return vj


def Density(mixture_ratio):
    # 推進剤の密度を求める
    p = (float)(1 + mixture_ratio) / \
            ((1 / H2_density) + (1 * mixture_ratio / O2_density))
    return p


def Finert(mixture_ratio):
    p = Density(mixture_ratio)
    finert_reverse = (1 / finerts - 1) * (p / ps) + 1
    finert = 1.0 / finert_reverse
    return finert


def Payload_Ratio(mixture_ratio,entrance_temp):
    vj = Vj(mixture_ratio,entrance_temp)
    finert = Finert(mixture_ratio)
    payload_ratio = (math.exp(-1.0 * delta_v / vj) - finert) / (1 - finert)
    return payload_ratio


def Thrust(mixture_ratio,entrance_temp):
    PR = Payload_Ratio(mixture_ratio,entrance_temp)
    return payload_mass * a_i / PR

#ガスジェネレーターサイクルのポンプ駆動動力
def gasg_pump_power(y,p):
    powp = p*pow(10,6)
    ans = (3.1 + 1.94*y)*1.2*powp/(1.14*1000*0.7) + (1 + 2.16*y)*1.6*20.64*pow(10,6)/(71*0.7)
    return ans

#ガスジェネレーターサイクルのタービン発生動力
def gasg_turbin_power(y,p):
    if p==0:
        return -1
    powp = p*pow(10,6)
    #ans = 4.1 * 8.09 * pow(10,3) * 0.5 * 780 * y * (1-pow(10000.0/powp,0.286))
    a = 1-pow(10000.0/powp,0.286)
    b = 4.1*8.09*1000*0.5*780
    return a*b*y


def cal_y(x):
    return x/(1-x)


def cal_xmin():
    xmin = 100000
    p_xmin = 0
    
    for x in range(0,10000):
        x = x*0.0001
        y = cal_y(x)
        for p in range(500,1000):
            p = p * 0.01
            if abs(gasg_pump_power(y,p) - gasg_turbin_power(y,p)) < 3 :
                if x < xmin:
                    xmin = x
                    p_xmin = p
                
    print("xmin:",xmin)
    print("p:",p_xmin)
    print("pump_power:",gasg_pump_power(cal_y(xmin),p_xmin))
    print("turbin_power:",gasg_turbin_power(cal_y(xmin),p_xmin))
    
if __name__ == '__main__':
    print("<前回課題>")
    tc = 300
    print("断熱火炎温度Tf:",Adiabatic_Temp(3.1, tc))
    print("排気速度:",Vj(3.1,tc))
    print("finert:",Finert(3.1))
    print("ペイロード比:",Payload_Ratio(3.1,tc))
    print("推力:",Thrust(3.1,tc))
    print("")
    
    
    print("<Gas Generator Cycle>")
    cal_xmin()
    gasg_tc = 300 #ガスジェエレーターサイクルの燃焼室燃焼前温度
    print("断熱火炎温度Tf:",Adiabatic_Temp(3.1, gasg_tc))
    print("排気速度:",Vj(3.1,tc))
    print("finert:",Finert(2.89))
    print("ペイロード比:",Payload_Ratio(2.89,gasg_tc))
    print("推力:",Thrust(2.89,gasg_tc))
    print("")
    
    print("<Staged Combustion Cycle>")
    stagec_tc = 838.64 #staged combustionの燃焼室燃焼前温度
    print("断熱火炎温度Tf:",Adiabatic_Temp(3.1, stagec_tc))
    print("排気速度:",Vj(3.1,stagec_tc))
    print("finert:",Finert(3.1))
    print("ペイロード比:",Payload_Ratio(3.1,stagec_tc))
    print("推力:",Thrust(3.1,stagec_tc))
    
    
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import math
mpl.rcParams['font.family'] = 'AppleGothic'
# 設定,要求したパラメータ
chamber_diameter = 0.25  # SSME 燃焼室直径 cd
chamber_pressure = 20.64 * math.pow(10, 6)  # N/m^2,SSME 燃焼室圧力 cp
tank_d = 2.5  # SSME燃料タンク直径 td
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


def mr_Change(mixture_ratio, mode):
    # 質量混合比をモルに変える moxider/mfuel
    if mode == 0:
        moll_ratio = mixture_ratio / 16
        # Oの原子量
    else:
        moll_ratio = mixture_ratio / 19
        # Fの原子量
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


def product_H2_1(moll_ratio):  # H2+F2=2HF
    if moll_ratio > 1:
        return 0
    else:
        return 1 - moll_ratio


def product_F2(moll_ratio):
    if moll_ratio > 1:
        return moll_ratio - 1
    else:
        return 0


def product_HF(moll_ratio):
    if moll_ratio > 1:
        return 2
    else:
        return 2 * moll_ratio


def Adiabatic_Temp(mixture_ratio, mode):
    if mode == 0:
        moll_ratio = mr_Change(mixture_ratio, 0)
    else:
        moll_ratio = mr_Change(mixture_ratio, 1)
    pH2_0 = product_H2_0(moll_ratio)
    pO2 = product_O2(moll_ratio)
    pH2o = product_H2o(moll_ratio)
    pH2_1 = product_H2_1(moll_ratio)
    pF2 = product_F2(moll_ratio)
    pHF = product_HF(moll_ratio)
    if mode == 0:  # oxider is O2
        enth = pH2_0 * cp + pO2 * cp + pH2o * cp
        hof = 241.9 * pH2o * 1000  # 単位J 生成熱
        temp = hof / enth + 300
        if temp > 2500:
            temp = (2 * hof + 726 * enth) / (enth + 0.0004 * hof)
        return temp
    else:  # oxider is F2
        enth = pH2_1 * cp + pF2 * cp + pHF * cp
        hof = 271.2 * pHF * 1000
        temp = (2 * hof + 300 * enth) / (enth + 0.0004 * hof)
        temp = hof / enth + 300
        if temp > 2500:
            temp = (2 * hof + 300 * enth) / (enth + 0.0004 * hof)
        return temp


def Products_mass(mixture_ratio, mode):
    # 水素1molとした時の反応後の重さ/mol
    moll_ratio = mr_Change(mixture_ratio, mode)
    pH2_0 = product_H2_0(moll_ratio)
    pO2 = product_O2(moll_ratio)
    pH2o = product_H2o(moll_ratio)
    pH2_1 = product_H2_1(moll_ratio)
    pF2 = product_F2(moll_ratio)
    pHF = product_HF(moll_ratio)
    if mode == 0:
        m = pH2_0 * 2 + pO2 * 32 + pH2o * 18 / (pH2_0 + pO2 + pH2o)
    else:
        m = pH2_1 * 2 + pF2 * 38 + pHF * 20 / (pH2_1 + pF2 + pHF)
    return m


def Vj(mixture_ratio, mode):
    m = Products_mass(mixture_ratio, mode)
    temp = Adiabatic_Temp(mixture_ratio, mode)
    gammaratio = (gamma - 1) / gamma
    powg = math.pow(PaP0ratio, gammaratio)
    vj = math.sqrt(2 * eta_j * (gamma / (gamma - 1))
                   * (R / m) * temp * (1 - powg))
    return vj


def Density(mixture_ratio, mode):
    # 推進剤の密度を求める
    if mode == 0:
        p = (float)(1 + mixture_ratio) / \
            ((1 / H2_density) + (1 * mixture_ratio / O2_density))
    else:
        p = (float)(1 + mixture_ratio) / \
            ((1 / H2_density) + (1 * mixture_ratio / F2_density))
    return p


def Finert(mixture_ratio, mode):
    p = Density(mixture_ratio, mode)
    finert_reverse = (1 / finerts - 1) * (p / ps) + 1
    finert = 1.0 / finert_reverse
    return finert


def Payload_Ratio(mixture_ratio, mode):
    vj = Vj(mixture_ratio, mode)
    finert = Finert(mixture_ratio, mode)
    payload_ratio = (math.exp(-1.0 * delta_v / vj) - finert) / (1 - finert)
    return payload_ratio


def Thrust(mixture_ratio, mode):
    PR = Payload_Ratio(mixture_ratio, mode)
    return payload_mass * a_i / PR


def Isp(mixture_ratio, mode):
    return Vj(mixture_ratio, mode) / g


def Thrust_Coefficient():  # Cf
    gammaratio_1 = (2.0 * math.pow(gamma, 2)) / (gamma - 1)
    gammaratio_2 = (gamma + 1) / (gamma - 1)
    gammaratio_3 = (gamma - 1) / gamma
    pow_1 = math.pow(2.0 / (gamma + 1), gammaratio_2)
    pow_2 = math.pow(PaP0ratio, gammaratio_3)
    CF2 = gammaratio_1 * pow_1 * (1 - pow_2)
    CF = math.sqrt(CF2)
    return CF


def Throat_A(mixture_ratio, mode):  # A*
    thrust = Thrust(mixture_ratio, mode)
    thrust_coefficient = Thrust_Coefficient()
    throat_a = thrust / (thrust_coefficient * chamber_pressure)
    return throat_a


def Chamber_Length(mixture_ratio, mode):
    throat_A = Throat_A(mixture_ratio, mode)
    adiabatic_temp = Adiabatic_Temp(mixture_ratio, mode)
    # 燃焼室断面積
    A1 = math.pi * math.pow(chamber_diameter, 2) / 4
    A1Atratio = A1 / throat_A  # 燃焼室断面積比の逆数
    #print("燃焼室断面積比：", A1Atratio)
    sqrt_1 = math.sqrt(Products_mass(mixture_ratio, mode) /
                       (gamma * R * adiabatic_temp))
    pow_1 = math.pow((gamma + 1) / 2, (gamma + 1) / 2 * (gamma - 1))
    L = (float)(1 / A1Atratio) * (ts / (sqrt_1 * pow_1))
    #print("燃焼室長さ:", L)
    return L


def Density_Liquid(mixture_ratio, mode):
    # Densityとorder違い
    if mode == 0:
        p = (1 + mixture_ratio) / ((1 / 70.8) + (1 * mixture_ratio / 1140))
    else:
        p = (1 + mixture_ratio) / ((1 / 70.8) + (1 * mixture_ratio / 1510))
    return p


def Tank_Length(mixture_ratio, mode):
    payload_ratio = Payload_Ratio(mixture_ratio, mode)
    finert = Finert(mixture_ratio, mode)
    density＿liquid = Density_Liquid(mixture_ratio, mode)
    propellant_mass = ((1 / payload_ratio) - 1) * (1 - finert) * payload_mass
    #print("propellant mass:", propellant_mass, "kg")
    tank_length = propellant_mass / \
        (density_liquid * (math.pi / 4) * tank_d * tank_d)
    return tank_length


def mrtemp(mode):
    x = np.arange(0.0, 15.0, 0.1)
    vAdiabatic_Temp = np.vectorize(Adiabatic_Temp)
    plt.plot(x, vAdiabatic_Temp(x, mode))
    plt.xlabel('mixture ratio')
    plt.ylabel('adiabatic flame temperature')
    plt.title('燃燒溫度')
    plt.savefig('燃焼温度.png')
    plt.show()


def mrpayload(mode):
    x = np.arange(0.0, 15.0, 0.1)
    vPayload_Ratio = np.vectorize(Payload_Ratio)
    plt.plot(x, vPayload_Ratio(x, mode))
    plt.xlabel('mixture ratio')
    plt.ylabel('payload ratio')
    plt.ylim(0, 0.1)
    plt.title('ペイロ-ド比')
    plt.savefig('ペイロード比.png')
    plt.show()


def mrisp(mode):
    x = np.arange(0.0, 15.0, 0.1)
    vIsp = np.vectorize(Isp)
    plt.plot(x, vIsp(x, mode))
    plt.xlabel('mixture ratio')
    plt.ylabel('Isp')
    plt.ylim(250, 475)
    plt.title('比推力')
    plt.savefig('比推力.png')
    plt.show()


def mrfinert(mode):
    x = np.arange(0.0, 15.0, 0.1)
    vFinert = np.vectorize(Finert)
    plt.plot(x, vFinert(x, mode))
    plt.xlabel('mixture ratio')
    plt.ylabel('Finert')
    plt.ylim(0.05, 0.35)
    plt.title('構造質量比')
    plt.savefig('構造質量比.png')
    plt.show()


def Max_Payload(mode):
    max_mr = 0
    max_payload = 0
    for i in range(1, 20):
        if Payload_Ratio(i, mode) > max_payload:
            max_mr = i
            max_payload = Payload_Ratio(i, mode)
        else:
            max_mr = max_mr
    for j in range((max_mr - 1) * 10, (max_mr + 1) * 10):
        if Payload_Ratio(1.0 * j / 10, mode) > max_payload:
            max_mr = 1.0 * j / 10
            max_payload = Payload_Ratio(1.0 * j / 10, mode)
        else:
            max_mr = max_mr
    adtmax = Adiabatic_Temp(max_mr, mode)
    ispmax = Isp(max_mr, mode)
    clmax = Chamber_Length(max_mr, mode)
    tlmax = Tank_Length(max_mr, mode)

    print("ペイロード比:", max_payload)
    print("燃焼室圧力:", chamber_pressure, "Pa")
    if mode == 0:
        print("推進剤:液体水素と液体酸素")
    else:
        print("推進剤:液体水素と液体フッ素")
    print("質量混合比:", max_mr)
    print("断熱火炎温度:", adtmax, "K")
    print("比推力:", ispmax)
    print("燃焼室径:", chamber_diameter, "m")
    print("燃焼室長さ:", clmax, "m")
    print("タンク径:", tank_d, "m")
    print('タンク長さ:', tlmax, "m")


if __name__ == '__main__':
    mode = input("choose O[enter:o] or F[enter:f]:")
    if mode == "o":  # 酸化剤が酸素
        #mrtemp(0)
        #mrpayload(0)
        #mrisp(0)
       # mrfinert(0)
        #Max_Payload(0)
        print("temp:",Adiabatic_Temp(3.3, 0))
        print("Vj:",Vj(3.3,0))
        print("Payload_ratio:",Payload_Ratio(3.3,0))
        print("Thrust:",Thrust(3.3,0))
    elif mode == "f":  # 酸化剤がフッ素
        mrtemp(1)
        mrpayload(1)
        mrisp(1)
        mrfinert(1)
        Max_Payload(1)
    else:
        print("ERROR")
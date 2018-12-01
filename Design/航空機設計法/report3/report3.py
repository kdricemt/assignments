#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 14 14:16:33 2018

@author: keidaiiiyama
"""

import math as m
import matplotlib.pyplot as plt

inch2feet = 0.0833
PI = 3.1415
rad2deg = 180 / PI
deg2rad = PI / 180
ft3togal = 7.480

#胴体計算
print("胴体計算")
#width
economy = 356
buisness = 64
subsec_ne = 5
subsec_nb = 3
seet_per_subsec_e = 6
seet_per_subsec_b = 4
seet_width_e = 21
seet_width_b = 30
corridor = 25
corridor_b = 30
wall = 5
subsec_we = seet_per_subsec_e * seet_width_e + corridor
subsec_wb = seet_per_subsec_b * seet_width_b + corridor_b
print("subsetion width(economy)",subsec_we)
print("subsetion width(buisness)",subsec_wb)
d_in1 = subsec_we * subsec_ne + wall * (subsec_ne - 1)
d_in2 = subsec_wb * subsec_nb + wall * (subsec_nb - 1)
d_out = d_in1 + 12
print("d_in(1st floor): ",d_in1,"[in]")
print("d_in(2nd floor): ",d_in2,"[in]")
print("d_out: ",d_out,"[in] ",d_out*0.0833,"[ft]")

#length
seetpitch_buisness = 65 #inch
seetpitch_economy = 33 #inch
seet_per_row_e = subsec_ne * seet_per_subsec_e
seet_per_row_b = subsec_nb * seet_per_subsec_b
row_economy = int(economy / seet_per_row_e + 1)
row_buisness = int(buisness / seet_per_row_b + 1)
print("seet_per_row(economy)",seet_per_row_e)
print("seet_per_row(buisness)",seet_per_row_b)
print("row_number (economy)",row_economy)
print("row_number (buisness)",row_buisness)
economy_length = seetpitch_economy * row_economy + corridor * 2
buisness_length = seetpitch_buisness * row_buisness + corridor * 2
print("economy_length",economy_length*0.0833,"[ft]")
print("buisness_length",buisness_length*0.0833,"[ft]")
lf = 1.3 * buisness_length * inch2feet #ft
print("length:",lf,"[ft]")
print()

# Thrust

T = 244044
print("推力計算")
print("Thrust:",T/3,"[Ib]")

#主翼計算
print()
print("主翼計算")
ar = 9.0
S = 7798.4
print("S",S,"[ft^2]")

b = m.sqrt(S*ar) #[ft]
print("wing span: ",b,"[ft]")
lamda = 0.3
w_to = 896000
w_pl = 93476
mff = 0.63
w_f = (1-mff)*w_to - w_pl
print("w_f:",w_f)
q = 194 #巡航時動圧
cla = 5 * deg2rad
alpha_i = 0 #胴体取り付け角
t_root = 0.18
t_tip = 0.12
rho_f = 6.7

#長さ
cr = 2/(1 + lamda) * m.sqrt(S/ar) #ft
print("C_root:",cr)
ct = 2/(1 + 1/lamda) * m.sqrt(S/ar) #ft
print("c_tip:",ct)
print()

#厚さ
thickness_root = 110 * t_root
thickness_tip = ct * t_tip
print("thickness_root",thickness_root,"[ft]")
print("thickness_tip",thickness_tip,"[ft]")
mode = "plot"

if mode=="plot":
 plt.rcParams["font.size"] = 18
 fig = plt.figure(figsize = (20,5))
 ax1 = fig.add_subplot(111)
 x = [-b,-80,-40,0,40,80,b]
 y1 = [thickness_tip,4.2,11.2,18,11.2,4.2,thickness_tip]
 ax1.set_xlabel("distance from root[ft]")
 ax1.set_ylabel("thickness[ft]")
 ln1 = ax1.plot(x,y1,color="blue",label="thickness")

 ax2 = ax1.twinx()
 y2 = [ct,30,70,100,70,30,ct]
 ax2.set_ylabel("chord[ft]")
 ln2 = ax2.plot(x,y2,color="red",label="chord")
 plt.legend()
 plt.savefig("./images/wingthickness.png")
 plt.show()
 print()

# 平均空力翼弦
c_mean = 2/3 * cr * (1+lamda+lamda*lamda)/(1+lamda)
print("cmean",c_mean)
cl_cruise = (w_to - 0.4 * w_f)/(q * S)
print("cl_cruise",cl_cruise)
alpha_i = cl_cruise/cla
print("alpha_i:",alpha_i)
v_t = 0.54 * S * S/b * t_root * (1 + lamda * m.sqrt(t_tip/t_root) + lamda **2 * t_tip/t_root) \
                                / ((1 + lamda)**2)
print("燃料タンク容量vt",v_t)
v_f = w_f/rho_f
print("必要な燃料体積vf",v_f)
print()

#wheel
print("タイヤのサイジング")
w_to = 896000
w_nw = 0.15 * w_to /4
print("w_nw: ",w_nw)
w_mw = 0.85 * w_to /8
print("w_mw: ",w_mw)
main_d = 1.63 * pow(w_mw,0.315)
print("main_d: ",main_d, "[in]  ", main_d * inch2feet, "[ft]")
main_width = 0.1043 * pow(w_mw,0.480)
print("main_width: ",main_width, "[in]  ", main_width * inch2feet, "[ft]")
nose_d = main_d * 0.7
print("nose_d: ",nose_d, "[in]  ", nose_d * inch2feet, "[ft]")
nose_width = main_width * 0.7
print("nose_width: ",nose_width, "[in]  ", nose_width * inch2feet, "[ft]")
print()

#Weight
print("重量の決定")
print("主翼")
n_ult = 3.8
r_angle = 30 #後退角
w_mzf = w_to - w_f
S_w = S * 1.2
print("S_w:", S_w)
print("w_mzf: ",w_mzf)
r_angle_half = m.atan((b / 2 * m.tan(deg2rad * r_angle) + ct/4 - cr/4)/(b/2))
print("r_angle_half: ",rad2deg * r_angle_half)
w_wing = 0.0017 * w_mzf * pow(b/m.cos(r_angle_half),0.75) * (1 + m.sqrt(6.25 / b * m.cos(r_angle_half))) \
       * pow(n_ult,0.55) * m.pow(b * S_w /(t_root * cr * w_mzf * m.cos(r_angle_half)),0.30)
print("w_wing: ",w_wing)

print("胴体")
s_fus = 2 * (12.2 * 65.5 + 45 * 12.2 + 65.5 * 45)
print("s_fus: ",s_fus," [ft^2]")
vd = 822 #[ft]
w_fus = 0.0065 * m.sqrt(vd) * 1.85 * pow(s_fus,1.2)
print("w_fus: ",w_fus," [Ib]")

print("ナセル重量")
t_to = 244044
w_nacelle = 0.065 * t_to
print("w_nacelle",w_nacelle)

print("脚重量")
k_gr = 1.0
w_mg = k_gr * (40 + 0.16*pow(w_to,0.75) + 0.019*w_to + 1.5 * pow(10,-5) * pow(w_to,1.5))
w_ng = k_gr * (20 + 0.10*pow(w_to,0.75) + 2.0 * pow(10,-6) * pow(w_to,1.5))
w_g = w_mg + w_ng;
print("w_mg: ",w_mg,"[Ib]")
print("w_ng: ",w_ng,"[Ib]")
print("w_g: ",w_g,"[Ib]")

print("推進系統重量")
w_eng = 15596 * 3
w_p = 1.16 * w_eng + 5950
print("w_eng: ",w_eng,"[Ib]")
print("w_p: ",w_p,"[Ib]")

print("装備品重量")
f_fix = 0.08
w_fix = f_fix * w_to
print("w_fix: ",w_fix,"[Ib]")

print("運用に必要なアイテム重量")
n_crew = 17
w_op = 187 * n_crew + 35 * 420
print("w_op: ",w_op,"[Ib]")

print("運用空虚重量")
w_tfo = w_to * 0.005
print("w_tfo: ",w_tfo,"[Ib]")
w_crew = 3485
print("w_crew: ",w_crew,"[Ib]")
w_me = w_wing + w_fus + w_nacelle + w_g
w_oe = w_me + w_p + w_fix + w_tfo + w_crew + w_op
print("w_oe: ",w_oe,"[Ib]")
w_oe2 = 471943 + w_tfo + w_crew + w_op
print("W_oe(report 2): ",w_oe2,"[Ib]")
print("error:", (w_oe2 - w_oe)/w_oe)

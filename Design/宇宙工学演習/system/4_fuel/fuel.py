#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import math

w_e = 0
w_en = 810 #タンクを含まない機器重量
w_t = 0 #タンク重量
w_d = 0 #ドライ重量
g = 9.8
deltav_sk = 399.06 #軌道姿勢制御
deltav_ap = 1869.1 #アポジ点でのキックモーター
isp_sk = 170
isp_ap = 280
sk_fuel = 0
ap_fuel = 0
i = 0

while(i<20):
    i += 1
    w_t = 0.1 * (sk_fuel + ap_fuel)
    w_e = w_en + w_t
    w_d = 1.07 * 1.17 * w_e
    sk_fuel = w_d * (math.exp(deltav_sk/(g * isp_sk))-1)
    ap_fuel = (w_d + sk_fuel) * (math.exp(deltav_ap/(g * isp_ap))-1)
    print("count:",i,"sk_fuel",sk_fuel,"ap_fuel",ap_fuel)
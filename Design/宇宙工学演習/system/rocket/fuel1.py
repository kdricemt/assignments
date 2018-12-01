import matplotlib.pyplot as plt
import math
import numpy as np

W_P2 = np.linspace(5, 200, 7000)
W_P1 = np.array([])


def cal_xi2(w_p2):
    return 0.975 * w_p2 / (3.059 + (1.145 - 0.02 * math.log10(w_p2)) * w_p2)


def cal_xi1(w_p1, w_p2):
    return 0.995 * w_p1 / (5.327 + (1.145 - 0.02 * math.log10(w_p2)) * w_p2 + (1.145 - 0.02 * math.log10(w_p1)) * w_p1)


def delta(w_p1, w_p2):
    return 11776 - 9.8 * (430 * math.log(1 / (1 - cal_xi1(w_p1, w_p2))) +
                          455 * math.log(1 / (1 - cal_xi2(w_p2))))


a_list = [10, 1000]

for j in W_P2:
    while(delta(a_list[0], j) >= 0 and delta(a_list[1], j) <= 0):
        if(abs(delta(a_list[1], j)) > 0.001):
            if(delta((a_list[0] + a_list[1]) / 2, j) > 0):
                a_list[0] = (a_list[0] + a_list[1]) / 2
            if(delta((a_list[0] + a_list[1]) / 2, j) < 0):
                a_list[1] = (a_list[0] + a_list[1]) / 2

        else:
            print(delta((a_list[0] + a_list[1]) / 2, j))
            print(a_list)
            W_P1 = np.append(W_P1, [(a_list[0] + a_list[1]) / 2])
            a_list = [10, 1000]
            break

plt.plot(W_P2, W_P1)
plt.xlabel("W_P2[t]")
plt.ylabel("W_P1[t]")
plt.show()

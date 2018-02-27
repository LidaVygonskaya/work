import numpy as np
from matplotlib import pyplot as plt
#initial data
#prostranstvo
x_0 = 0
x_N = 500 #meters
N = 100
h = (x_N - x_0) / N
V_i = h #odnomernayua i ravnomernaya

#plast
fi = 0.2 #poristost
k = (9.868233 * (10 ** (-13))) * 10 ** (-3)
atm = 101325.0
P_init = 80.0 * atm #Pa
P_0 = 130 * atm
P_Nx = 30 * atm
ro_0 = 1000 #kg/m3 new
P_01 = 80 * atm #new
P_02 = 80 * atm # new
fi_0 = 0.2

#fluid
c_r = (10 ** (-5)) / atm # new
c_f = (10 ** (-4)) / atm  # szhimaemost
mu = 10 ** (-3)  # Pa * s

#stepping
tau = 86400.0 #s- time step
time_max = 365.25 * 4 * tau
dt_out_fields = 10 * tau # frecuency of writing


#setka po x
x = np.arange(x_0, x_N, (x_N - x_0) / N)
x = np.append(x, x_N)
x = list(x)


#coefficients

def count_a(index, x_list):
    T_X_plus = (k / mu) * (1 / (x_list[index + 1] - x_list[index]))
    T_X_minus = (k / mu) * (1 / (x_list[index] - x_list[index - 1]))
    return -(T_X_plus + T_X_minus + (h * fi * c_f / tau))


def count_b(index, x_list):
    return (k / mu) * (1 / (x_list[index + 1] - x_list[index]))


def count_c(index, x_list):
    return (k / mu) * (1 / (x_list[index] - x_list[index - 1]))


def count_d(P_i_n):
    return -(P_i_n * h * fi * c_f / tau)


def thomas_method(list_a, list_b, list_c, list_d):
    w_1 = list_a[0]
    q_list = []
    w_list = [w_1]
    for j in range(1, len(list_a)):
        q = list_b[j - 1] / w_list[j - 1]
        q_list.append(q)
        w = list_a[j] - list_c[j] * q_list[j - 1]
        w_list.append(w)

    g_list = [list_d[0] / w_list[0]]
    for j in range(1, len(list_d)):
        g = (list_d[j] - list_c[j] * g_list[j - 1]) / w_list[j]
        g_list.append(g)

    list_P = [g_list[-1]]
    for j in range(len(g_list) - 2, -1, -1):
        P = g_list[j] - q_list[j] * list_P[-1]
        list_P.append(P)
    # Add usloviya na granitsah
    list_P.append(P_0)
    list_P.reverse()
    list_P.append(P_Nx)
    return list_P


def initialize_P(P0, Pn, Pinit):
    list_P = []
    list_P.append(P0)
    for i in range(1, N):
        list_P.append(Pinit)
    list_P.append(Pn)
    return list_P


def initialize_d(list_P, list_c, list_b):
    list_d = []
    for j in range(1, N):
        d = count_d(list_P[j])
        list_d.append(d)

    list_d[0] -= list_c[0] * list_P[0]
    list_d[-1] -= list_b[-1] * list_P[-1]
    return list_d


#zavisimost plotnosti ot davleniya
def ro(p):
    return ro_0 * (1.0 + c_f * (p - P_02))


#zavisimost fi ot davleniya
def fi_p(p):
    return fi_0 * (1.0 + c_r * (p - P_01))


def beta(p):
    return fi_0 * c_r * ro(p) + c_f * ro_0 * fi_p(p)


def count_c_nonlinear(list_P, index, x_list):
    return (k / mu) * (1 / (x_list[index] - x_list[index - 1])) * ro(max(list_P[index], list_P[index - 1]))


def count_b_nonlinear(list_P, index, x_list):
    return (k / mu) * (1 / (x_list[index + 1] - x_list[index])) * ro(max(list_P[index], list_P[index + 1]))


def count_a_nonlinear(list_P, index, x_list):
    T_X_plus = (k / mu) * (1 / (x_list[index + 1] - x_list[index])) * ro(max(list_P[index], list_P[index + 1]))
    T_X_minus = (k / mu) * (1 / (x_list[index] - x_list[index - 1])) * ro(max(list_P[index], list_P[index - 1]))
    return -(T_X_minus + T_X_plus + h * beta(list_P[index]) / tau)


def count_d_nonlinear(P_i_n):
    return -(P_i_n * h * beta(P_i_n) / tau)


def initialize_d_nonlinear(list_P, list_c, list_b):
    list_d = []
    for j in range(1, N):
        d = count_d_nonlinear(list_P[j])
        list_d.append(d)

    list_d[0] -= list_c[0] * list_P[0]
    list_d[-1] -= list_b[-1] * list_P[-1]
    return list_d


P_list = initialize_P(P_0, P_Nx, P_init)
a_list = []
b_list = []
c_list = []

for i in range(1, N):
    a = count_a_nonlinear(P_list, i, x)
    b = count_b_nonlinear(P_list, i, x)
    c = count_c_nonlinear(P_list, i, x)
    a_list.append(a)
    b_list.append(b)
    c_list.append(c)

d_list = initialize_d_nonlinear(P_list, c_list, b_list)

time = tau
counter = 1

while time <= time_max:
    P_list = thomas_method(a_list, b_list, c_list, d_list)
    if counter % 10 == 0:
        name = 'graph_nonlinear ' + str(counter)
        plt.plot(x, P_list)
        plt.xlabel('x')
        plt.ylabel('P')
        plt.savefig(name)
        plt.clf()
    counter += 1
    time += tau
    a_list = []
    b_list = []
    c_list = []
    for i in range(1, N):
        a = count_a_nonlinear(P_list, i, x)
        b = count_b_nonlinear(P_list, i, x)
        c = count_c_nonlinear(P_list, i, x)
        a_list.append(a)
        b_list.append(b)
        c_list.append(c)
    d_list = initialize_d_nonlinear(P_list, c_list, b_list)




# #LINEAR PART
# a_list = []
# b_list = []
# c_list = []
#
# for i in range(1, N):
#     a = count_a(i, x)
#     b = count_b(i, x)
#     c = count_c(i, x)
#
#     a_list.append(a)
#     b_list.append(b)
#     c_list.append(c)
#
# P_list = initialize_P(P_0, P_Nx, P_init)
# d_list = initialize_d(P_list, c_list, b_list)
#
#
# time = tau
# counter = 1
#
# while time <= time_max:
#     P_list = thomas_method(a_list, b_list, c_list, d_list)
#     if counter % 10 == 0:
#         name = 'graph ' + str(counter)
#         plt.plot(x, P_list)
#         plt.xlabel('x')
#         plt.ylabel('P')
#         plt.savefig(name)
#         plt.clf()
#     counter += 1
#     time += tau
#     d_list = initialize_d(P_list, c_list, b_list)




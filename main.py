import numpy as np
import math
from matplotlib import pyplot as plt
import copy
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

#Newtons method
delta_0 = 1000.0 #nachalnoe priblizhenie
delta_max = 10 ** (-3) #zadannaya tochnost
delta_gr = 1000.0

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


def thomas_method(list_a, list_b, list_c, list_d, P0, PN):
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
    list_P.append(P0)
    list_P.reverse()
    list_P.append(PN)
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
    t_x_plus = (k / mu) * (1 / (x_list[index + 1] - x_list[index])) * ro(max(list_P[index], list_P[index + 1]))
    t_x_minus = (k / mu) * (1 / (x_list[index] - x_list[index - 1])) * ro(max(list_P[index], list_P[index - 1]))
    return -(t_x_minus + t_x_plus + h * beta(list_P[index]) / tau)


def count_d_nonlinear(P_i_n):
    return -(P_i_n * h * beta(P_i_n) / tau)


def initialize_d_nonlinear(list_P, list_c, list_b):
    list_d = []
    for j in range(1, N):
        d = count_d_nonlinear(list_P[j])
        list_d.append(d)

    #list_d[0] -= list_c[0] * list_P[0]
    #list_d[-1] -= list_b[-1] * list_P[-1]
    return list_d


def initialize_abc(list_P, x_list):
    list_a = []
    list_b = []
    list_c = []
    for j in range(1, N):
        _a = count_a_nonlinear(list_P, j, x_list)
        _b = count_b_nonlinear(list_P, j, x_list)
        _c = count_c_nonlinear(list_P, j, x_list)
        list_a.append(_a)
        list_b.append(_b)
        list_c.append(_c)
    list_abc = [list_a, list_b, list_c]
    return list_abc


#Newtons method
def count_discrepancy(list_abcd, list_p, index):
    a, b, c, d = list_abcd
    return c[index] * list_p[index] + a[index] * list_p[index + 1]  + b[index] * list_p[index + 2] - d[index]


def diff_simple(list_a, list_b, list_c):
    list_diff = []
    for j in range(N - 1):
        list_diff_i = [list_c[j], list_a[j], list_b[j]]
        list_diff.append(list_diff_i)
    return list_diff


def diff(list_p_k, list_delta, list_abcd):
    list_diff = []
    for j in range(N - 1):
        R_minus = list_abcd[2][j] * (list_p_k[j] + list_delta[j]) + list_abcd[0][j] * list_p_k[j + 1] + list_abcd[1][j] * list_p_k[j + 2]

        R_0 = list_abcd[2][j] * list_p_k[j] + list_abcd[0][j] * (list_p_k[j + 1] + list_delta[j + 1]) + list_abcd[1][j] * list_p_k[j + 2]

        R_plus = list_abcd[2][j] * list_p_k[j] + list_abcd[0][j] * list_p_k[j + 1] + list_abcd[1][j] * (list_p_k[j + 2] + list_delta[j + 2])

        R_diff_minus = (R_minus - count_discrepancy(list_abcd, list_p_k, j)) / list_delta[j]
        R_diff_0 = (R_0 - count_discrepancy(list_abcd, list_p_k, j)) / list_delta[j + 1]
        R_diff_plus = (R_plus - count_discrepancy(list_abcd, list_p_k, j)) / list_delta[j + 2]
        list_diff_j = [R_diff_minus, R_diff_0, R_diff_plus]
        list_diff.append(list_diff_j)
    return list_diff


def check_diag(a, b, c, index):
    if math.fabs(a) <= math.fabs(b) + math.fabs(c):
        print("ne shoditsa " + str(index))


def count_norm(list_delta):
    list_delta_abs = [math.fabs(elem) for elem in list_delta[1:-1]]
    return max(list_delta_abs)


#NEWTONS METHOD
P_list_k = initialize_P(P_0, P_Nx, P_init) #init
P_list_n = initialize_P(P_0, P_Nx, P_init)
time = tau
counter = 1

while time < time_max:
    delta_list_k = [delta_gr for i in range(len(x))]
    for i in range(1, len(delta_list_k)):
        delta_list_k[i] = delta_0

    while count_norm(delta_list_k) > delta_max:
        print("norm " + str(count_norm(delta_list_k)))
        a_list = []
        b_list = []
        c_list = []
        for i in range(1, N):
            a = count_a_nonlinear(P_list_k, i, x)
            b = count_b_nonlinear(P_list_k, i, x)
            c = count_c_nonlinear(P_list_k, i, x)
            check_diag(a, b, c, i)

            a_list.append(a)
            b_list.append(b)
            c_list.append(c)

        d_list = initialize_d_nonlinear(P_list_n, c_list, b_list)
        abcd_list = [a_list, b_list, c_list, d_list]

        R_list = []
        for i in range(len(a_list)):
            R = count_discrepancy(abcd_list, P_list_k, i)
            R_list.append(R)

        #R_diff = diff(P_list_k, delta_list_k, abcd_list)
        R_diff = diff_simple(a_list, b_list, c_list)
        a_delta_list = []
        b_delta_list = []
        c_delta_list = []

        for i in range(N - 1):
            a_delta = R_diff[i][1]
            b_delta = R_diff[i][2]
            c_delta = R_diff[i][0]
            check_diag(a_delta, b_delta, c_delta, i)

            a_delta_list.append(a_delta)
            b_delta_list.append(b_delta)
            c_delta_list.append(c_delta)

        c_delta_list[0] = 0.0
        b_delta_list[-1] = 0.0

        d_delta_list = [- elem for elem in R_list]

        delta_list_k = thomas_method(a_delta_list, b_delta_list, c_delta_list, d_delta_list, delta_gr, delta_gr)
        for i in range(1, len(P_list_k) - 1):
            P_list_k[i] = P_list_k[i] + delta_list_k[i]

    P_list_n = copy.deepcopy(P_list_k)
    if counter % 10 == 0:
        name = 'graph_' + str(counter)
        plt.plot(x, P_list_n)
        plt.xlabel('x')
        plt.ylabel('y')
        plt.savefig(name)
        plt.clf()
    counter += 1
    time += tau


#NONLINEAR PART
# P_list = initialize_P(P_0, P_Nx, P_init)
# a_list = []
# b_list = []
# c_list = []
#
# for i in range(1, N):
#     a = count_a_nonlinear(P_list, i, x)
#     b = count_b_nonlinear(P_list, i, x)
#     c = count_c_nonlinear(P_list, i, x)
#     a_list.append(a)
#     b_list.append(b)
#     c_list.append(c)
#
# d_list = initialize_d_nonlinear(P_list, c_list, b_list)
#
# time = tau
# counter = 1
#
# while time <= time_max:
#     P_list = thomas_method(a_list, b_list, c_list, d_list)
#     if counter % 10 == 0:
#         name = 'graph_nonlinear ' + str(counter)
#         plt.plot(x, P_list)
#         plt.xlabel('x')
#         plt.ylabel('P')
#         plt.savefig(name)
#         plt.clf()
#     counter += 1
#     time += tau
#     a_list = []
#     b_list = []
#     c_list = []
#     for i in range(1, N):
#         a = count_a_nonlinear(P_list, i, x)
#         b = count_b_nonlinear(P_list, i, x)
#         c = count_c_nonlinear(P_list, i, x)
#         a_list.append(a)
#         b_list.append(b)
#         c_list.append(c)
#     d_list = initialize_d_nonlinear(P_list, c_list, b_list)




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




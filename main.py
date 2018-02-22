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

#fluid
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


a_list = []
b_list = []
c_list = []

for i in range(1, N):
    a = count_a(i, x)
    b = count_b(i, x)
    c = count_c(i, x)

    a_list.append(a)
    b_list.append(b)
    c_list.append(c)


#srazu posle t = 0
d_list = []
for i in range(1, N):
    d = count_d(P_init)
    d_list.append(d)

d_list[0] = d_list[0] - c_list[0] * P_0
d_list[-1] = d_list[-1] - b_list[-1] * P_Nx

#thomas method
w_1 = a_list[0]
q_list = []
w_list = [w_1]
for i in range(1, len(a_list)):
    q = b_list[i - 1] / w_list[i - 1]
    q_list.append(q)
    w_i = a_list[i] - c_list[i] * q_list[i - 1]
    w_list.append(w_i)

g_list = [d_list[0] / w_list[0]]
for i in range(1, len(d_list)):
    g = (d_list[i] - c_list[i] * g_list[i - 1]) / w_list[i]
    g_list.append(g)

print(len(g_list))

P_list = [g_list[-1]]
for i in range(len(g_list) - 2, -1, -1):
    P = g_list[i] - q_list[i] * P_list[-1]
    P_list.append(P)



#Add usloviya na granitsah
P_list.append(P_0)
P_list.reverse()
P_list.append(P_Nx)


time = tau

while time <= time_max:
    
    time += tau






plt.plot(x, P_list)
plt.xlabel('x')
plt.ylabel('P')
plt.show()

#plt.savefig('name.png')


#
#
#
#
#
#
# for time in range(int(time_max / tau)):
#     d_list = []
#     for i in range(1, N):
#         if time == 0:
#             d = count_d(P_init)
#
#         d_list.append(d)
#
#     d_list = [] # from i  = 1 to i = N - 1
#     d = cou
#
#



import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from matplotlib import pyplot as plt
import copy
#initial data
#prostranstvo
x_0 = 0.0
x_N = 500.0 #meters
y_0 = 0.0
y_N = 500.0
N_x = 100
N_y = 100
N = 100.0
#h_x = (x_N - x_0) / (N_x - 1)
#h_y = (y_N - y_0) / (N_y - 1)
h = (x_N - x_0) / (N_x - 1) #v sluchae esli x i y gorizontalniye
#h = (x_N - x_0) / N
V_ij = h ** 3.0 #odnomernayua i ravnomernaya


#plast
fi = 0.2 #poristost
k = (9.868233 * (10 ** (-13))) * 10 ** (-3)

k_x = k
k_y = k
atm = 101325.0
P_init = 80.0 * atm #Pa
P_0 = 130 * atm
P_Nx = 30 * atm
ro_0 = 1000 #kg/m3 new
P_01 = 80 * atm #new
P_02 = 80 * atm # new
fi_0 = 0.2
z = h
#fluid
c_r = (10 ** (-5)) / atm # new
c_f = (10 ** (-4)) / atm  # szhimaemost
mu = 10 ** (-3)  # Pa * s

#skvazhina
r_w = 0.108 # meters
P_w = 30 * atm # Pa
x_s = 250
y_s = 250


#stepping
tau = 86400.0 #s- time step
#time_max = 365.25 * 4 * tau
time_max = 365.25 * 1 * tau
dt_out_fields = 10 * tau # frecuency of writing


#setka po x
x = np.arange(x_0, x_N, (x_N - x_0) / N)
x = np.append(x, x_N)
x = list(x)

delta_0 = 1000.0 #nachalnoe priblizhenie
delta_max = 10 ** (-3) #zadannaya tochnost
delta_gr = 1000.0


def ro(p):
    return ro_0 * (1.0 + c_f * (p - P_02))


def fi_p(p):
    return fi_0 * (1.0 + c_r * (p - P_01))


def count_re():
    return h * math.exp(-math.pi / 2.0)


def count_tolerance(matrix):
    matrix_abs = np.absolute(matrix)
    return np.amax(matrix_abs)

#dP/dx = 0 na granitsah
def count_tx(p):
    tx_plus = np.zeros((N_x - 1, N_y - 1))
    tx_minus = np.zeros((N_x - 1, N_y - 1))
    for i in range(N_y - 1):
        for j in range(N_x - 1):
            if j == N_x - 2:
                tx_plus[i][j] = 0.0
            else:
                tx_plus[i][j] = (k_x * (h ** 2) / (mu * h)) * ro(max(p[i][j], p[i][j + 1]))

        for j in range(N_x - 1):
            if j == 0:
                tx_minus[i][j] = 0.0
            else:
                tx_minus[i][j] = (k_x * (h ** 2) / (mu * h)) * ro(max(p[i][j], p[i][j - 1]))
    tx_list = [tx_minus, tx_plus]
    return tx_list


def count_ty(p):
    ty_plus = np.zeros((N_x - 1, N_y - 1))
    ty_minus = np.zeros((N_x - 1, N_y - 1))
    for j in range(N_x - 1):
        for i in range(N_y - 1):
            if i == N_y - 2:
                ty_plus[i][j] = 0.0
            else:
                ty_plus[i][j] = (k_y * (h ** 2) / (mu * h)) * ro(max(p[i][j], p[i + 1][j]))
        for i in range(N_y - 1):
            if i == 0:
                ty_minus[i][j] = 0.0
            else:
                ty_minus[i][j] = (k_y * (h ** 2) / (mu * h)) * ro(max(p[i][j], p[i - 1][j]))
    ty_list = [ty_minus, ty_plus]
    return ty_list



#P = const na granitsah
# def count_tx(p_gr):
#     tx_plus = np.zeros((N_x - 1, N_y - 1))
#     tx_minus = np.zeros((N_x - 1, N_y - 1))
#     for i in range(N_y - 1):
#         for j in range(N_x - 1):
#             tx_plus[i][j] = (k_x * (h ** 2) / (mu * h)) * ro(max(p_gr[i + 1][j + 1], p_gr[i + 1][j + 2]))
#         for j in range(N_x - 1):
#             tx_minus[i][j] = (k_x * (h ** 2) / (mu * h)) * ro(max(p_gr[i + 1][j + 1], p_gr[i + 1][j]))
#     tx_list = [tx_minus, tx_plus]
#     return tx_list
#
#
# def count_ty(p_gr):
#     ty_plus = np.zeros((N_x - 1, N_y - 1))
#     ty_minus = np.zeros((N_x - 1, N_y - 1))
#     for j in range(N_x - 1):
#         for i in range(N_y - 1):
#             ty_plus[i][j] = (k_y * (h ** 2) / (mu * h)) * ro(max(p_gr[i + 1][j + 1], p_gr[i + 2][j + 1]))
#         for i in range(N_y - 1):
#             ty_minus[i][j] = (k_y * (h ** 2) / (mu * h)) * ro(max(p_gr[i + 1][j + 1], p_gr[i][j + 1]))
#     ty_list = [ty_minus, ty_plus]
#     return ty_list


def count_beta(p):
    beta = np.zeros((N_y - 1, N_x - 1))
    for i in range(N_y - 1):
        for j in range(N_x - 1):
            beta[i][j] = fi_0 * c_r * ro(p[i][j]) + c_f * ro_0 * fi_p(p[i][j])
    return beta


def count_wi(p):
    wi = np.zeros((N_x - 1, N_y - 1))
    re = count_re()
    for i in range(N_y - 1):
        for j in range(N_x - 1):
            if i == 45 and j == 45:
                wi[i][j] = -(ro(p[i][j]) * 2 * math.pi * k * h / mu) * (1 / math.log(re / r_w))
    return wi


def count_r(p_gr, c, g, a, f, b, d):
    r = np.zeros((N_x - 1, N_y - 1))
    for i in range(N_y - 1):
        r[i] = c[i] * p_gr[i + 1][:-2] + g[i] * p_gr[i][1:-1] + a[i] * p_gr[i + 1][1:-1] + f[i] * p_gr[i + 2][1:-1] + b[i] * p_gr[i + 1][2:] - d[i] # Check!!!!!?????
    return r


def check_diag(a, b, c, g, f):
    for i in range(N_y - 1):
        for j in range(N_x - 1):
            if math.fabs(a[i][j]) < math.fabs(b[i][j]) + math.fabs(c[i][j]) + math.fabs(g[i][j]) + math.fabs(f[i][j]):
                print("ne shod " + str(i) + str(j))


p = np.ones((N_x - 1, N_y - 1)) * P_init
p_gr = np.ones((N_y + 1, N_x + 1)) * P_init
delta = np.ones((N_y + 1, N_x + 1)) * delta_0
delta_k = np.zeros((N_y + 1, N_x + 1))
ksi = delta.copy()
ksi_k = np.ones((N_y + 1, N_x + 1))
q_list = []
time_list = []


time = tau
counter = 1
while time < time_max:
    delta = np.ones((N_y + 1, N_x + 1)) * delta_0
    delta_k = np.zeros((N_y + 1, N_x + 1))
    p_n = p.copy()
    while count_tolerance(delta_k - delta) > delta_max:
        ksi = np.ones((N_y + 1, N_x + 1)) * delta_0
        beta = count_beta(p)
        wi = count_wi(p)
        c, b = count_tx(p)#!!!!!!!!!!!!!!!!!!!!!
        g, f = count_ty(p)#!!!!!!!!!!!!!!!!!
        a = -(b + c + g + f) + wi - V_ij * beta / tau
        d = wi * P_w - p_n * (V_ij * beta / tau)

        check_diag(a, b, c, g, f)

        r_x_diff_minus = c.copy()
        r_x_diff_plus = b.copy()
        r_y_diff_minus = g.copy()
        r_y_diff_plus = f.copy()
        r_diff = a.copy()
        r = count_r(p_gr, c, g, a, f, b, d)
        delta_k = delta.copy()

        while count_tolerance(ksi - ksi_k) > delta_max:
            ksi_k = ksi.copy()
            for i in range(N_y - 1):
                ksi[i + 1][1:-1] = (1 / r_diff[i]) * (
                -r[i] - r_x_diff_plus[i] * ksi[i + 1][2:] - r_y_diff_plus[i] * ksi[i + 2][1:-1] - r_x_diff_minus[i] *
                ksi[i + 1][:-2] - r_y_diff_minus[i] * ksi[i][1:-1])
            tolerance = count_tolerance(ksi - ksi_k)
            #print(tolerance)

        delta = ksi.copy()
        tolerance = count_tolerance(delta - delta_k)
        print(str(tolerance) + " meow")


        for i in range(N_y - 1):
            for j in range(N_x - 1):
                p[i][j] += delta[i + 1][j + 1]


        p_gr += delta






    if counter % 1 == 0:
        name = 'graph_' + str(counter)
        fig = plt.figure()
        ax = fig.gca(projection='3d')

        # Make data.
        X = np.arange(x_0 + h / 2, x_N, h)
        Y = np.arange(y_0 + h / 2, y_N, h)
        X, Y = np.meshgrid(X, Y)

        # Plot the surface.
        surf = ax.plot_surface(X, Y, p, cmap=cm.coolwarm,
                               linewidth=0, antialiased=False)
        plt.savefig(name)
        plt.clf()

    q = wi[45][45] * (p[45][45] - P_w)
    q_list.append(q)
    time_list.append(time)
    print(q_list)
    print(time_list)

    time += tau
    counter += 1

# file = open('q_time.txt', 'w')
# for i in range(len(q_list)):
#     file.write(str(time_list[i]) + '  ' + str(q_list[i]) + '\n')
# file.close()
#
# plt.plot(time_list, q_list)
# plt.show()
#
#
#



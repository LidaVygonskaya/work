import math
import matplotlib.pyplot as plt
import numpy as np


class Layer:
    def __init__(self):
        self.atm = 101325.0
        self.ro_oil_0 = 900.0
        self.ro_water_0 = 1000.0
        self.P_01 = 80 * self.atm  # new
        self.P_02 = 80 * self.atm
        self.c_f_oil = (10.0 ** (-4)) / self.atm
        self.c_f_water = (10.0 ** (-4)) / self.atm
        self.c_r = (10.0 ** (-5)) / self.atm
        self.fi_0 = 0.2
        self.mu_oil = 10.0 * (10.0 ** (-3))
        self.mu_water = 10.0 ** (-3)
        self.k = (9.868233 * (10 ** (-13))) * 10 ** (-3)
        self.s_water_init = 0.2
        self.s_oil_init = 1 - self.s_water_init
        self.pressure_cap_init = [] #list
        self.pressure_oil_init = 80.0 * self.atm
        self.pressure_water_init = []
        for m in range(len(self.pressure_cap_init[1: -1])):
            self.pressure_water_init.append(self.pressure_oil_init - self.pressure_cap_init[1: -1][m])



        #left side
        self.pressure_water_left = 130.0 * self.atm
        self.s_water_left = 1.0
        self.s_oil_left = 0.0
        #self.pressure_oil_left = self.pressure_water_left + self.pressure_cap_init[0]

        #right side
        self.pressure_water_right = self.pressure_water_init
        self.pressure_oil_right = self.pressure_oil_init
        self.s_oil_right = self.s_oil_init
        self.s_water_right = self.s_water_init

        #space
        self.x_0 = 0.0
        self.x_N = 500.0  # meters
        self.N = 100
        self.h = (self.x_N - self.x_0) / (self.N - 1)
        self.z = 10.0
        self.V_ij = self.h ** 2.0 * self.z

    @staticmethod
    def count_k_r(_s_water):
        k_r_water = _s_water ** 2.0
        k_r_oil = (1.0 - _s_water) ** 2.0
        return [k_r_water, k_r_oil]

    def count_ro_water_oil(self,_pressure_water, _pressure_oil):
        ro_water = self.ro_water_0 * (1.0 + self.c_f_water * (_pressure_water - self.P_02))
        ro_oil = self.ro_oil_0 * (1.0 + self.c_f_oil* (_pressure_oil - self.P_02))
        return [ro_water, ro_oil]

    def count_ro_water(self, _pressure_water):
        return self.ro_water_0 * (1.0 + self.c_f_water * (_pressure_water - self.P_02))

    def count_ro_oil(self, _pressure_oil):
        return self.ro_oil_0 * (1.0 + self.c_f_oil* (_pressure_oil - self.P_02))

    def count_fi(self, pressure_oil):
        return self.fi_0 * (1.0 + self.c_r * (pressure_oil - self.P_01))

    @staticmethod
    def count_pcap(p_cap_graph, s_water_list):
        p_cap_list = []
        s_w_graph = list(p_cap_graph.keys())
        for i in range(len(s_water_list)):
            for j in range(len(s_w_graph)):
                if s_water_list[i] < s_w_graph[j]:
                    p_cap = p_cap_graph.get(s_w_graph[j - 1]) + (
                                p_cap_graph.get(s_w_graph[j]) - p_cap_graph.get(
                            s_w_graph[j - 1])) / (s_w_graph[j] - s_w_graph[j - 1]) * (
                                        s_water_list[i] - s_w_graph[j - 1])
                    p_cap_list.append(p_cap)
                    break
        return p_cap_list




class Solver:
    def __init__(self):
        self.delta_0 = 1000.0
        self.delta_max = 10 ** (-3)
        self.tau = 8640.0  # s- time step
        # time_max = 365.25 * 4 * tau
        self.time_max = 365.25 * 1 * self.tau

    @staticmethod
    def count_t_water_minus(layer,  pressure_water, s_water, index):
        ro_water = layer.count_ro_water(max(pressure_water[index - 1], pressure_water[index]))
        k_r_water = layer.count_k_r(max(s_water[index - 1], s_water[index]))[0]
        t_x_water_minus = (layer.k * k_r_water / layer.mu_water) * (1 / layer.h) * ro_water
        return t_x_water_minus

    @staticmethod
    def count_t_water_plus(layer, pressure_water, s_water, index):
        ro_water = layer.count_ro_water(max(pressure_water[index], pressure_water[index + 1]))
        k_r_water = layer.count_k_r(max(s_water[index], s_water[index + 1]))[0]
        t_x_water_plus = (layer.k * k_r_water / layer.mu_water) * (1 / layer.h) * ro_water
        return t_x_water_plus

    @staticmethod
    def count_t_oil_minus(layer, pressure_oil, s_water, index):
        ro_oil = layer.count_ro_oil(max(pressure_oil[index - 1], pressure_oil[index]))
        k_r_oil = layer.count_k_r(max(s_water[index - 1], s_water[index]))[1]
        t_x_oil_minus = (layer.k * k_r_oil / layer.mu_oil) * (1 / layer.h) * ro_oil
        return t_x_oil_minus

    @staticmethod
    def count_t_oil_plus(layer, pressure_oil, s_water, index):
        ro_oil = layer.count_ro_oil(max(pressure_oil[index], pressure_oil[index + 1]))
        k_r_oil = layer.count_k_r(max(s_water[index], s_water[index + 1]))[1]
        t_x_oil_plus = (layer.k * k_r_oil / layer.mu_oil) * (1 / layer.h) * ro_oil
        return t_x_oil_plus

    @staticmethod
    def count_c1_p(layer, fi_n, s_water, ro_water, index ):
        return (fi_n[index] * s_water[index] * layer.ro_water_0 * layer.c_f_water) + (s_water[index] * ro_water[index] * layer.c_r * layer.fi_0)

    @staticmethod
    def count_c2_p(layer, fi_n, s_water, ro_oil, index):
        return (1.0 - s_water[index]) * (fi_n[index] * layer.ro_oil_0 * layer.c_f_oil + layer.c_r * layer.fi_0 * ro_oil[index])


    @staticmethod
    def count_b(ro_water, ro_oil, t_x_water_plus, t_x_oil_plus, index):
        A = ro_oil[index] / ro_water[index]
        b = A * t_x_water_plus + t_x_oil_plus
        return b

    @staticmethod
    def count_c(ro_oil, ro_water, t_x_water_minus, t_x_oil_minus, index):
        A = ro_oil[index] / ro_water[index]
        c = A * t_x_water_minus + t_x_oil_minus
        return c

    @staticmethod
    def count_a(ro_water, ro_oil, C1p, C2p, b, c, index):
        A = ro_oil[index] / ro_water[index]
        return -(b + c + A * C1p + C2p)

    @staticmethod
    def count_d(pressure_oil,pressure_oil_n, pressure_cap, ro_water, ro_oil, t_x_water_minus, t_x_water_plus,a, b, c,C1p, C2p, index):
        A = ro_oil[index] / ro_water[index]
        r = b[index] * pressure_oil[index + 2] + c[index] * pressure_oil[index] + a[index] * pressure_oil[index + 1] - A * t_x_water_plus * (pressure_cap[index + 2] - pressure_cap[index + 1]) - A * t_x_water_minus * (pressure_cap[index] - pressure_cap[index + 1]) + (A * C1p + C2p) * pressure_oil_n[index + 1]
        return -r

    @staticmethod
    def count_norm(list_delta):
        list_delta_abs = [math.fabs(elem) for elem in list_delta[1:-1]]
        return max(list_delta_abs)

    @staticmethod
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


layer = Layer()
solver = Solver()
pressure_cap_graph = {} #from graph {s_w: pressure_cap}
file = open('Pcap(Sw).txt', 'r')
for line in file.readlines():
    line = line.rstrip()
    s = line.split('\t')
    pressure_cap_graph.update({float(s[0]): float(s[1]) * layer.atm})



#initialize S_water
s_water = []
for i in range(layer.N - 2):
    s_water.append(layer.s_water_init)

#initialize S_oil
s_oil = []
for i in range(layer.N - 2):
    s_oil.append(1 - layer.s_water_init)

#initialize Oil pressure without gr
P_oil = []
for j in range(layer.N - 2):
    P_oil.append(layer.pressure_oil_init)

#initialize Pressure_cap from 1 to N - 1
P_cap = []
s_water_graph = list(pressure_cap_graph.keys())
for i in range(len(s_water)):
    for j in range(len(s_water_graph)):
        if s_water[i] < s_water_graph[j]:
            p_cap = pressure_cap_graph.get(s_water_graph[j - 1]) + (pressure_cap_graph.get(s_water_graph[j]) - pressure_cap_graph.get(s_water_graph[j - 1])) / (s_water_graph[j] - s_water_graph[j - 1]) * (s_water[i] - s_water_graph[j-1])
            P_cap.append(p_cap)
            break

#Water pressure
P_water = []
for i in range (len(P_oil)):
    P_water.append(P_oil[i] - P_cap[i])

#Adding gr conditions
#s_water
s_water = [layer.s_water_left] + s_water
s_water.append(layer.s_water_right)

#s_oil
s_oil = [1 - layer.s_water_left] + s_oil
s_oil.append(1 - layer.s_water_right)

#P_cap
P_cap = [pressure_cap_graph.get(s_water[0])] + P_cap
P_cap.append(pressure_cap_graph.get(s_water[-1]))
print(P_cap)

#P_water
P_water = [layer.pressure_water_left] + P_water
P_water.append(P_water[-1])


#P_oil
P_oil = [P_water[0] + P_cap[0]] + P_oil
P_oil.append(P_oil[-1])

#fi_list_n

time = solver.tau
counter = 1
s_water_from_oil = []
while time < solver.time_max:

    delta_list_k = [solver.delta_0 for i in range(layer.N)]
    P_oil_n = P_oil.copy()
    s_water_n = s_water.copy()
    s_oil_n = s_oil.copy()
    fi_list_n = [layer.count_fi(P_oil_n[i]) for i in range(1, layer.N - 1)]
    P_cap = [pressure_cap_graph.get(s_water_n[0])] + layer.count_pcap(pressure_cap_graph, s_water_n[1:-1])
    P_cap.append(pressure_cap_graph.get(s_water_n[-1]))
    print(len(P_cap))


    while solver.count_norm(delta_list_k) > solver.delta_max:

        #print(solver.count_norm(delta_list_k))
        a_list = []
        b_list = []
        c_list = []
        d_list = []
        t_water_minus_list = []
        t_water_plus_list = []
        t_oil_minus_list = []
        t_oil_plus_list = []
        ro_water = []
        ro_oil = []
        fi_list = []

        for i in range(1, layer.N - 1):
            t_water_minus_list.append(solver.count_t_water_minus(layer, P_water, s_water, i))
            t_oil_minus_list.append(solver.count_t_oil_minus(layer, P_oil, s_water, i))
            t_water_plus_list.append(solver.count_t_water_plus(layer, P_water, s_water, i))
            t_oil_plus_list.append(solver.count_t_oil_plus(layer, P_oil, s_water, i))
            ro_water.append(layer.count_ro_water(P_water[i]))
            ro_oil.append(layer.count_ro_oil(P_oil[i]))
            fi_list.append(layer.count_fi(P_oil[i]))

        for i in range(layer.N - 2):
            b_list.append(solver.count_b(ro_water, ro_oil, t_water_plus_list[i], t_oil_plus_list[i], i))
            c_list.append(solver.count_c(ro_water, ro_oil, t_water_minus_list[i], t_oil_minus_list[i], i))
            c1_p = solver.count_c1_p(layer, fi_list_n, s_water[1:-1], ro_water, i)
            c2_p = solver.count_c2_p(layer, fi_list_n, s_water[1:-1], ro_oil, i)
            a_list.append(solver.count_a(ro_water, ro_oil, c1_p, c2_p, b_list[i], c_list[i], i))
            d_list.append(solver.count_d(P_oil,P_oil_n, P_cap, ro_water, ro_oil, t_water_minus_list[i], t_water_plus_list[i], a_list, b_list, c_list, c1_p, c2_p, i))

        delta_list_k = solver.thomas_method(a_list, b_list, c_list, d_list, solver.delta_0, solver.delta_0)
        for i in range(1, layer.N - 1):
            P_oil[i] = P_oil[i] + delta_list_k[i]

    for i in range(layer.N - 2):
        s_water[i + 1] = s_water_n[i + 1] + (solver.tau/(fi_list[i] * ro_water[i])) * (t_water_plus_list[i] * (P_oil[i + 2] - P_oil_n[i + 1]) + t_water_minus_list[i] * (P_oil[i] - P_oil[i + 1]))
        s_oil[i + 1] = s_oil_n[i + 1] + ((solver.tau/(fi_list[i] * ro_oil[i]))) * (t_oil_plus_list[i] * (P_oil[i + 2] - P_oil[i + 1] + P_cap[i + 2] - P_cap[i + 1]) + t_oil_minus_list[i] * (P_oil[i] - P_oil[i + 1] + P_cap[i] - P_cap[i + 1]))

    s_water_from_oil = [1 - s_oil[i] for i in range(len(s_oil))]
    print(solver.count_norm([s_water[i] - s_water_n[i] for i in range(len(s_water))]))

    if solver.count_norm([s_water[i] - s_water_n[i] for i in range(len(s_water))]) > 0.1:
        solver.tau = solver.tau / 2.0
        P_oil = P_oil_n.copy()
        s_water = s_water_n.copy()
        fi_list = fi_list_n.copy()

    else:
        x = np.arange(layer.x_0, layer.x_N, layer.h)
        x = np.append(x, layer.N)
        k_r_water = [s_water[i] ** 2 for i in range(len(s_water))]
        k_r_oil = [(1 - s_water[i]) ** 2 for i in range(len(s_water))]
        plt.plot(x, P_oil)
        plt.show()
        plt.plot(s_water, k_r_water)
        plt.plot(s_water, k_r_oil)
        plt.show()

        time += solver.tau
        counter += 1
    #Count S n + 1

























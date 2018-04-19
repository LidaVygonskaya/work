import math
import matplotlib.pyplot as plt
import numpy as np

from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Layer import Layer
from oil_water_filtration.Solver import SolverSlau


class Solver:
    def __init__(self):
        self.delta_0 = 1000.0
        self.delta_max = 10 ** (-3)
        self.tau_default = 86400.0  # s- time step
        self.tau = self.tau_default
        # time_max = 365.25 * 4 * tau
        self.time_max = 365.25 * 1 * 86400 #self.tau

    @staticmethod
    def count_t_water_minus(layer,  pressure_water, s_water, index):
        ro_wat = layer.count_ro_water(max(pressure_water[index - 1], pressure_water[index]))
        k_r_water = layer.count_k_r(max(s_water[index - 1], s_water[index]))[0]
        t_x_water_minus = (layer.k * k_r_water / layer.mu_water) * (1 / layer.h) ** 2.0 * ro_wat
        return t_x_water_minus

    @staticmethod
    def count_t_water_plus(layer, pressure_water, s_water, index):
        ro_wat = layer.count_ro_water(max(pressure_water[index], pressure_water[index + 1]))
        k_r_water = layer.count_k_r(max(s_water[index], s_water[index + 1]))[0]
        t_x_water_plus = (layer.k * k_r_water / layer.mu_water) * (1 / layer.h) ** 2.0 * ro_wat
        return t_x_water_plus

    @staticmethod
    def count_t_oil_minus(layer, pressure_oil, s_water, index):
        _ro_oil = layer.count_ro_oil(max(pressure_oil[index - 1], pressure_oil[index]))
        k_r_oil = layer.count_k_r(max(s_water[index - 1], s_water[index]))[1]
        t_x_oil_minus = (layer.k * k_r_oil / layer.mu_oil) * (1 / layer.h) ** 2.0 * _ro_oil
        return t_x_oil_minus

    @staticmethod
    def count_t_oil_plus(layer, pressure_oil, s_water, index):
        _ro_oil = layer.count_ro_oil(max(pressure_oil[index], pressure_oil[index + 1]))
        k_r_oil = layer.count_k_r(max(s_water[index], s_water[index + 1]))[1]
        t_x_oil_plus = (layer.k * k_r_oil / layer.mu_oil) * (1 / layer.h) ** 2.0 * _ro_oil
        return t_x_oil_plus

    @staticmethod
    def count_c1_p(layer, solver, fi_n, s_water_n, ro_water, index):
        return ((fi_n[index] * s_water_n[index] * layer.ro_water_0 * layer.c_f_water) + (s_water_n[index] * ro_water[index] * layer.c_r * layer.fi_0)) / solver.tau

    @staticmethod
    def count_c2_p(layer, solver,  fi_n, s_water_n, ro_oil, index):
        return ((1.0 - s_water_n[index]) * (fi_n[index] * layer.ro_oil_0 * layer.c_f_oil + layer.c_r * layer.fi_0 * ro_oil[index])) / solver.tau


    @staticmethod
    def count_b(solverSlau, ro_water, ro_oil, t_x_water_plus, t_x_oil_plus, pressure_oil, index):
        A = ro_oil[index] / ro_water[index]
        #A = 1.0
        b = A * t_x_water_plus + t_x_oil_plus
        solverSlau.set_matrix_coefficients(index, index + 1, pressure_oil, b)
        return b

    @staticmethod
    def count_c(solverSlau, ro_oil, ro_water, t_x_water_minus, t_x_oil_minus, pressure_oil,  index):
        A = ro_oil[index] / ro_water[index]
        #A = 1.0
        c = A * t_x_water_minus + t_x_oil_minus

        solverSlau.set_matrix_coefficients(index, index - 1, pressure_oil, c)

        # solverSlau.add_coefficient(index, index, -c)
        # solverSlau.add_nevyaz(index, -c * pressure_oil[index])
        # solverSlau.add_nevyaz(index, c * pressure_oil[index + 1])
        # if index != 0:
        #     solverSlau.add_coefficient(index, index - 1, c)
        return c

    @staticmethod
    def count_a(solverSlau, ro_water, ro_oil, C1p, C2p, pressure_oil, pressure_oil_n,  index):
        A = ro_oil[index] / ro_water[index]
        #A = 1.0
        #a = -(b + c + A * C1p + C2p)
        solverSlau.set_matrix_coeffitients_a_d(index, index, pressure_oil, pressure_oil_n, A * C1p + C2p)
        #solverSlau.add_coefficient(index, index, -A * C1p + C2p)
        #solverSlau.add_nevyaz(index, A * C1p + C2p * pressure_oil[index + 1])
        return 0

    @staticmethod
    def count_d(solverSlau, pressure_oil_n, pressure_cap, ro_water, ro_oil, t_x_water_minus, t_x_water_plus, C1p, C2p, index):
        A = ro_oil[index] / ro_water[index]
        #A = 1.0
        #r = b[index] * pressure_oil[index + 2] + c[index] * pressure_oil[index] + a[index] * pressure_oil[index + 1] - A * t_x_water_plus * (pressure_cap[index + 2] - pressure_cap[index + 1]) - A * t_x_water_minus * (pressure_cap[index] - pressure_cap[index + 1]) + (A * C1p + C2p) * pressure_oil_n[index + 1]
        r_ost = - A * t_x_water_plus * (pressure_cap[index + 2] - pressure_cap[index + 1]) - A * t_x_water_minus * (pressure_cap[index] - pressure_cap[index + 1])# + (A * C1p + C2p) * pressure_oil_n[index + 1]
        solverSlau.add_nevyaz(index, -r_ost)
        return 0

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

cell_container = CellContainer(layer.N, layer)
for cell in cell_container.get_cells():
    index = cell_container.get_cells().index(cell)
    if index == 0:
        cell.cell_states[0].set_pressure_oil(layer.)
    elif index == 





pressure_cap_graph = {} #from graph {s_w: pressure_cap}
file = open('Pcap(Sw).txt', 'r')
for line in file.readlines():
    line = line.rstrip()
    s = line.split('\t')
    pressure_cap_graph.update({float(s[0]): float(s[1]) * layer.atm})



#initialize S_water
s_water = []
#s_water1 = np.matrix((1, layer.N))

for i in range(layer.N - 2):
    s_water.append(layer.s_water_init)
    #s_water1 = np.append(s_water1, layer.s_water_init)
#initialize S_oil
s_oil = []
#s_oil1 = np.matrix((1, layer.N))
for i in range(layer.N - 2):
    s_oil.append(1 - layer.s_water_init)
    #np.append(s_oil1, 1 - layer.s_water_init)

#initialize Oil pressure without gr
P_oil = []
#P_oil1 = np.matrix((1, layer.N))
for j in range(layer.N - 2):
    P_oil.append(layer.pressure_oil_init)

#initialize Pressure_cap from 1 to N - 1
P_cap = []
s_water_graph = list(pressure_cap_graph.keys())
for i in range(len(s_water)):
    for j in range(len(s_water_graph)):
        if s_water[i] <= s_water_graph[j]:
            p_cap = pressure_cap_graph.get(s_water_graph[j - 1]) + ((pressure_cap_graph.get(s_water_graph[j]) - pressure_cap_graph.get(s_water_graph[j - 1])) / (s_water_graph[j] - s_water_graph[j - 1])) * (s_water[i] - s_water_graph[j - 1])
            P_cap.append(p_cap)
            break

#Water pressure
P_water = []
for i in range(len(P_oil)):
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
print(len(P_cap))

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

solver_slau = SolverSlau(layer.N - 2)

while time < solver.time_max:

    delta_list_k = [0.0] + [solver.delta_0 for i in range(layer.N - 2)] + [0.0]
    P_oil_n = P_oil.copy()
    P_water_n = P_water.copy()
    s_water_n = s_water.copy()
    s_oil_n = s_oil.copy()
    fi_list_n = [layer.count_fi((P_oil_n[i] + P_water_n[i]) / 2.0) for i in range(1, layer.N - 1)]                     # 0.5 * (P_w+P_o)
    ro_water_n = [layer.count_ro_water((P_oil_n[i] + P_water_n[i]) / 2.0) for i in range(1, layer.N - 1)]
    P_cap = [pressure_cap_graph.get(s_water_n[0])] + layer.count_pcap(pressure_cap_graph, s_water_n[1:])
    #P_cap.append(pressure_cap_graph.get(s_water_n[-1]))
    print(len(P_cap))

    while solver.count_norm(delta_list_k) > solver.delta_max:
        #print(solver.count_norm(delta_list_k))
        # if solver.count_norm([delta_list_k[i] / P_oil_n[i] for i in range(len(delta_list_k))]) > 0.1:
        #     print("P_n+1 - P_n > 0.1")
        #     solver.tau = solver.tau / 2.0
        #     P_oil = P_oil_n.copy()
        #     fi_list = fi_list_n.copy()
        #     delta_list_k = [0.0] + [solver.delta_0 for i in range(layer.N - 2)] + [0.0]

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
            #b_list.append(solver.count_b(solver_slau, ro_water, ro_oil, t_water_plus_list[i], t_oil_plus_list[i], i))
            #c_list.append(solver.count_c(solver_slau, ro_oil, ro_water, t_water_minus_list[i], t_oil_minus_list[i], i))
            c1_p = solver.count_c1_p(layer, solver, fi_list_n, s_water_n[1:-1], ro_water, i)
            c2_p = solver.count_c2_p(layer, solver, fi_list_n, s_water_n[1:-1], ro_oil, i)
            #a_list.append(solver.count_a(solver_slau, ro_water, ro_oil, c1_p, c2_p, b_list[i], c_list[i], i))
            #d_list.append(solver.count_d(solver_slau, P_oil, P_oil_n, P_cap, ro_water, ro_oil, t_water_minus_list[i], t_water_plus_list[i], a_list, b_list, c_list, c1_p, c2_p, i))


        #count coefficients with solverSlau

        for i in range(solver_slau.e_count):
            solver.count_b(solver_slau, ro_water, ro_oil, t_water_plus_list[i], t_oil_plus_list[i],P_oil, i)
            solver.count_c(solver_slau, ro_oil, ro_water, t_water_minus_list[i], t_oil_minus_list[i], P_oil, i)
            solver.count_a(solver_slau, ro_water, ro_oil, c1_p, c2_p, P_oil, P_oil_n, i)
            solver.count_d(solver_slau, P_oil_n, P_cap, ro_water, ro_oil, t_water_minus_list[i], t_water_plus_list[i], c1_p, c2_p, i)




            #solver_slau.nevyaz_vector[i] = d_list[i]

        #delta_list_k = solver.thomas_method(a_list, b_list, c_list, d_list, 0.0, 0.0)
        solver_slau.solve_thomas_method(0.0, 0.0)
        delta_list_k_new = solver_slau.get_result()

        for i in range(1, layer.N - 1):
            P_oil[i] = P_oil[i] + delta_list_k[i]
            P_water[i] = P_oil[i] - P_cap[i]

    # Recalculate water density and transmissibilities
    t_water_minus_list = []
    t_water_plus_list = []
    for i in range(1, layer.N - 1):
        t_water_minus_list.append(solver.count_t_water_minus(layer, P_water, s_water, i))
        t_water_plus_list.append(solver.count_t_water_plus(layer, P_water, s_water, i))

    for i in range(layer.N - 2):
        ro_water[i] = layer.count_ro_water(P_water[i + 1])
        fi_list[i] = layer.count_fi(P_oil[i + 1])
        s_water[i + 1] = s_water_n[i + 1]*( fi_list_n[i] * ro_water_n[i])/(fi_list[i] * ro_water[i]) + (solver.tau/(fi_list[i] * ro_water[i])) * (t_water_plus_list[i] * (P_oil[i + 2] - P_oil[i + 1] + P_cap[i + 1] - P_cap[i + 2]) + t_water_minus_list[i] * (P_oil[i] - P_oil[i + 1] + P_cap[i + 1] - P_cap[i])) # + t_water_plus_list[i] * (P_cap[i + 1] - P_cap[i + 2]) + t_water_minus_list[i] * (P_cap[i + 1] - P_cap[i]))
        #s_oil[i + 1] = s_oil_n[i + 1] + (solver.tau/(fi_list[i] * ro_oil[i])) * (t_oil_plus_list[i] * (P_oil[i + 2] - P_oil[i + 1] + P_cap[i + 2] - P_cap[i + 1]) + t_oil_minus_list[i] * (P_oil[i] - P_oil[i + 1] + P_cap[i] - P_cap[i + 1]))

    s_water[layer.N - 1] = s_water[layer.N - 2]
    #s_water_from_oil = [1 - s_oil[i] for i in range(len(s_oil))]
    print(solver.count_norm([s_water[i] - s_water_n[i] for i in range(len(s_water))]))

    if solver.count_norm([s_water[i] - s_water_n[i] for i in range(len(s_water))]) > 0.01:
        print("S_n+1 - S_n > 0.1")
        solver.tau = solver.tau / 2.0
        P_oil = P_oil_n.copy()
        P_water = P_water_n.copy()
        s_water = s_water_n.copy()
        s_oil = s_oil_n.copy()
        ro_water = ro_water_n.copy()
        fi_list = fi_list_n.copy()

    else:

        #plt.plot(s_water, k_r_water)
        #plt.plot(s_water, k_r_oil)

        time += solver.tau
        solver.tau = max(solver.tau * 2.0, solver.tau_default)
        if 125 * 86400.0 - time - solver.tau <= 0:
            solver.tau = 125 * 86400.0 - time #- solver.tau

        if (abs(time - 86400.0*125) < 1.e-16):
            x = np.arange(layer.x_0, layer.x_N, layer.h)
            x = np.append(x, layer.x_N)
            #
            # plt.plot(x, P_water, 'bo', x, P_water, 'k')
            # plt.title(str(time / 86400.0) + 'days')
            # plt.xlabel('x')
            # plt.ylabel('P_water')
            # plt.show()
            # plt.plot(x, s_water, 'bo', x, s_water, 'k')
            # plt.title(str(time / 86400.0) + 'days')
            # plt.xlabel('x')
            # plt.ylabel('S_water')
            # plt.show()
            file = open('s(x)_count_125_days.txt', 'w')
            for i in range(len(x)):
                file.write(str(x[i]) + ' ' + str(s_water[i]) + '\n')
            file.close()

        counter += 1
    #Count S n + 1

























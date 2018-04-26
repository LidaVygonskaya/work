import math
import matplotlib.pyplot as plt
import numpy as np

from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
from oil_water_filtration.Solver import SolverSlau
from oil_water_filtration.OilWater import OilWater

class Solver:
    def __init__(self):
        self.delta_0 = 1000.0
        self.delta_max = 10 ** (-3)
        self.tau_default = 86400.0  # s- time step
        self.tau = self.tau_default
        # time_max = 365.25 * 4 * tau
        self.time_max = 365.25 * 1 * 86400 #self.tau

    @staticmethod
    #del
    def count_t_water_minus(layer,  pressure_water, s_water, index):
        ro_wat = layer.count_ro_water(max(pressure_water[index - 1], pressure_water[index]))
        k_r_water = layer.count_k_r(max(s_water[index - 1], s_water[index]))[0]
        t_x_water_minus = (layer.k * k_r_water / layer.mu_water) * (1 / layer.h) ** 2.0 * ro_wat
        return t_x_water_minus

    @staticmethod
    #del
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

    #Считаем во всех клетках кроме граничных


    @staticmethod
    def count_c2_p(layer, solver,  fi_n, s_water_n, ro_oil, index):
        return ((1.0 - s_water_n[index]) * (fi_n[index] * layer.ro_oil_0 * layer.c_f_oil + layer.c_r * layer.fi_0 * ro_oil[index])) / solver.tau


    @staticmethod
    def count_b(solverSlau, cell_flow):
        left_cell = cell_flow.get_left_cell()
        right_cell = cell_flow.get_right_cell()
        #_cell = cell_flow.get_left_cell()
        state_n_plus = left_cell.get_cell_state_n_plus()
        A = state_n_plus.get_ro_oil() / state_n_plus.get_ro_water()
        #b = A * t_x_water_plus + t_x_oil_plus
        b = A * cell_flow.t_oil_water[1] + cell_flow.t_oil_water[0]
        solverSlau.set_matrix_coefficients(cell_flow.get_left_cell().get_eq_index()[0], cell_flow.get_right_cell().get_eq_index()[0], left_cell.get_cell_state_n_plus().get_pressure_oil(), right_cell.get_cell_state_n_plus().get_pressure_oil(), b)
        return b

    @staticmethod
    def count_c(solverSlau, cell_flow):
        left_cell = cell_flow.get_left_cell()
        right_cell = cell_flow.get_right_cell()
        #_cell = cell_flow.get_right_cell()
        state_n_plus = right_cell.get_cell_state_n_plus()
        A = state_n_plus.get_ro_oil() / state_n_plus.get_ro_water()
        #c = A * t_x_water_minus + t_x_oil_minus
        c = A * cell_flow.t_oil_water[1] + cell_flow.t_oil_water[0]
        r_ost = - A * cell_flow.t_oil_water[1] * (left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        solverSlau.add_nevyaz(index, -r_ost)
        #solverSlau.set_matrix_coefficients(index, index - 1, state_n_plus.get_pressure_oil(), c)
        solverSlau.set_matrix_coefficients(cell_flow.get_right_cell().get_eq_index()[0], cell_flow.get_left_cell().get_eq_index()[0], right_cell.get_cell_state_n_plus().get_pressure_oil(), left_cell.get_cell_state_n_plus().get_pressure_oil(), c)
        return c

    @staticmethod
    def count_a(solverSlau, cell):
        #A = ro_oil[index] / ro_water[index]
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        A = state_n_plus.get_ro_oil() / state_n_plus.get_ro_water()
        solverSlau.set_matrix_coeffitients_a_d(cell.get_eq_index()[0], cell.get_eq_index()[0], state_n_plus.get_pressure_oil(), state_n.get_pressure_oil(), A * state_n_plus.get_c1_p() + state_n_plus.get_c2_p())


    @staticmethod
    def count_d(solverSlau, pressure_cap, t_x_water_minus, t_x_water_plus, cell, index):
        state_n_plus = cell.get_cell_state_n_plus()
        A = state_n_plus.get_ro_oil() / state_n_plus.get_ro_water()
        r_ost = - A * t_x_water_plus * (pressure_cap[index + 2] - pressure_cap[index + 1]) - A * t_x_water_minus * (pressure_cap[index] - pressure_cap[index + 1])
        solverSlau.add_nevyaz(index, -r_ost)
        return 0

    @staticmethod
    def count_norm(list_delta):
        list_delta_abs = [math.fabs(elem) for elem in list_delta[1:-1]]
        return max(list_delta_abs)


layer = Layer()
solver = Solver()
oil_water = OilWater()
cell_container = CellContainer(layer.N, layer)
cell_container.initialize_cells()

flow_array = []
for i in range(cell_container.get_len() - 1):
    cells = cell_container.get_cells()
    flow_array.append(Flow(cells[i], cells[i + 1]))



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

    #same
    P_oil_n = P_oil.copy()
    P_water_n = P_water.copy()
    s_water_n = s_water.copy()
    s_oil_n = s_oil.copy()
    fi_list_n = [layer.count_fi((P_oil_n[i] + P_water_n[i]) / 2.0) for i in range(1, layer.N - 1)]
    ro_water_n = [layer.count_ro_water((P_oil_n[i] + P_water_n[i]) / 2.0) for i in range(1, layer.N - 1)]
    #same

    P_cap = [pressure_cap_graph.get(s_water_n[0])] + layer.count_pcap(pressure_cap_graph, s_water_n[1:])

    #new
    for cell in cell_container.get_cells():
        cell.get_cell_state_n().set_equals_to(cell.get_cell_state_n_plus())
    #new

    while solver.count_norm(delta_list_k) > solver.delta_max:
        print (solver.count_norm(delta_list_k))
        solver_slau.set_zero()

        t_water_minus_list = []
        t_water_plus_list = []
        t_oil_minus_list = []
        t_oil_plus_list = []
        ro_water = []
        ro_oil = []
        fi_list = []

        for i in range(cell_container.get_len() - 1):
            cells = cell_container.get_cells()
            cell = cells[i]
            cell_state_n_plus = cell.get_cell_state_n_plus()
            for component in cell.layer.components:
                component_index = cell.layer.components.index(component)
                ro = component.count_ro(cell_state_n_plus.get_components_pressure()[component_index])
                k_r = component.count_k_r(cell_state_n_plus.get_components_saturation()[component_index])
                cell_state_n_plus.set_ro(ro, component_index)
                cell_state_n_plus.set_k_r(k_r, component_index)
            fi = cell.layer.count_fi(cell_state_n_plus.get_pressure_oil())
            cell_state_n_plus.set_fi(fi)

        for flow in flow_array:
            flow.count_cells_flows()

        for i in range(1, layer.N - 1):
            t_water_minus_list.append(solver.count_t_water_minus(layer, P_water, s_water, i))
            t_oil_minus_list.append(solver.count_t_oil_minus(layer, P_oil, s_water, i))
            t_water_plus_list.append(solver.count_t_water_plus(layer, P_water, s_water, i))
            t_oil_plus_list.append(solver.count_t_oil_plus(layer, P_oil, s_water, i))
            ro_water.append(layer.count_ro_water(P_water[i]))
            ro_oil.append(layer.count_ro_oil(P_oil[i]))
            fi_list.append(layer.count_fi(P_oil[i]))


        for flow in flow_array:
            solver.count_b(solver_slau, flow)
            solver.count_c(solver_slau, flow)

        for cell in cell_container.get_cells():
            cell.get_cell_state_n_plus().set_c1_p(cell.layer.count_c1_p_new(solver, cell))
            cell.get_cell_state_n_plus().set_c2_p(cell.layer.count_c2_p_new(solver, cell))
            solver.count_a(solver_slau, cell)


        for i in range(solver_slau.e_count):
            cell = cell_container.get_cells()[i + 1]
            solver.count_d(solver_slau, P_cap, t_water_minus_list[i], t_water_plus_list[i], cell, i)


        solver_slau.solve_thomas_method(0.0, 0.0)
        delta_list_k = solver_slau.get_result()
        solver_slau.clear_result()


        #adding delta and changing pressure
        for i in range(1, layer.N - 1):
            #same
            P_oil[i] = P_oil[i] + delta_list_k[i]
            P_water[i] = P_oil[i] - P_cap[i]

            #new the same
            cell = cell_container.get_cells()[i]
            state = cell.get_cell_state_n_plus()
            state.set_pressure_oil(state.get_pressure_oil() + delta_list_k[i])
            state.set_pressure_water(state.get_pressure_oil() - P_cap[i])


    # Recalculate water density and transmissibilities
    #transmissibilities same
    t_water_minus_list = []
    t_water_plus_list = []
    for i in range(1, layer.N - 1):
        t_water_minus_list.append(solver.count_t_water_minus(layer, P_water, s_water, i))
        t_water_plus_list.append(solver.count_t_water_plus(layer, P_water, s_water, i))


    #recount_density and fi water
    for i in range(1, cell_container.get_len() - 1):
        cells = cell_container.get_cells()
        cell = cells[i]
        cell_state = cell.get_cell_state_n_plus()
        for component in cell.layer.components:
            component_index = cell.layer.components.index(component)
            ro = component.count_ro(cell_state_n_plus.get_components_pressure()[component_index])
            k_r = component.count_k_r(cell_state_n_plus.get_components_saturation()[component_index])
            cell_state_n_plus.set_ro(ro, component_index)
            cell_state_n_plus.set_k_r(k_r, component_index)
        fi = cell.layer.count_fi(cell_state.get_pressure_oil())


    # recount transmissibilities new the same
    t_new_water = []
    for i in range(cell_container.get_len() - 1):
        cells = cell_container.get_cells()
        flow = Flow(cells[i], cells[i + 1])
        t_new_water.append(flow.count_cells_flows()[1])

    t_water_minus_list_new = [0.0] + t_new_water[:-1] + [0.0]
    t_water_plus_list_new = [0.0] + t_new_water[1:] + [0.0]


    #recount density fi and s
    for i in range(layer.N - 2):
        ro_water[i] = layer.count_ro_water(P_water[i + 1])
        fi_list[i] = layer.count_fi(P_oil[i + 1])
        s_water[i + 1] = s_water_n[i + 1] * (fi_list_n[i] * ro_water_n[i])/(fi_list[i] * ro_water[i]) + (solver.tau/(fi_list[i] * ro_water[i])) * (t_water_plus_list[i] * (P_oil[i + 2] - P_oil[i + 1] + P_cap[i + 1] - P_cap[i + 2]) + t_water_minus_list[i] * (P_oil[i] - P_oil[i + 1] + P_cap[i + 1] - P_cap[i]))
    s_water[layer.N - 1] = s_water[layer.N - 2]

    for i in range(1, cell_container.get_len() - 1):
        cells = cell_container.get_cells()
        cell_i = cells[i]
        cell_i_plus = cells[i + 1]
        cell_i_minus = cells[i - 1]

        cell_state_n_plus_i = cell_i.get_cell_state_n_plus()
        cell_state_n_plus_i_plus = cell_i_plus.get_cell_state_n_plus()
        cell_state_n_plus_i_minus = cell_i_minus.get_cell_state_n_plus()
        cell_state_n_i = cell_i.get_cell_state_n()

        coeff1 = cell_state_n_i.get_fi() * cell_state_n_i.get_ro_water() / (cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_water())
        coeff2 = t_water_plus_list_new[i] * (cell_state_n_plus_i_plus.get_pressure_oil() - cell_state_n_plus_i.get_pressure_oil() + P_cap[i] - P_cap[i + 1]) + t_water_minus_list_new[i] * (cell_state_n_plus_i_minus.get_pressure_oil() - cell_state_n_plus_i.get_pressure_oil() + P_cap[i] - P_cap[i - 1])#!!!!!!!!!!!!!!!!!!!!!!add to t_water
        s_water_new = cell_state_n_i.get_s_water() * coeff1 + (solver.tau / (cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_water())) * coeff2
        cell_state_n_plus_i.set_s_water(s_water_new)
        cell_state_n_plus_i.set_s_oil(1.0 - s_water_new)

    #cells = cell_container.get_cells()
    #cells[-2].get_state_n = cells[-1].

    if solver.count_norm([cell.get_cell_state_n_plus().get_s_water() - cell.get_cell_state_n().get_s_water() for cell in cell_container.get_cells()]) > 0.01:
    #if solver.count_norm([s_water[i] - s_water_n[i] for i in range(len(s_water))]) > 0.01:
        print("S_n+1 - S_n > 0.1")
        solver.tau = solver.tau / 2.0
        #same
        P_oil = P_oil_n.copy()
        P_water = P_water_n.copy()
        s_water = s_water_n.copy()
        s_oil = s_oil_n.copy()
        ro_water = ro_water_n.copy()
        fi_list = fi_list_n.copy()

        #the same
        for cell in cell_container.get_cells():
            cell.get_cell_state_n_plus().set_equals_to(cell.get_cell_state_n())

    else:
        time += solver.tau
        solver.tau = max(solver.tau * 2.0, solver.tau_default)
        if 125 * 86400.0 - time - solver.tau <= 0:
            solver.tau = 125 * 86400.0 - time #- solver.tau

        if abs(time - 86400.0 * 125) < 1.e-16:
            x = np.arange(layer.x_0, layer.x_N, layer.h)
            x = np.append(x, layer.x_N)
            #
            plt.plot(x, P_water, 'bo', x, P_water, 'k')
            plt.title(str(time / 86400.0) + 'days')
            plt.xlabel('x')
            plt.ylabel('P_water')
            plt.show()
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


























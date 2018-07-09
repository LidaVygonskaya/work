import numpy as np
import math
from oil_water_filtration.OilWater import OilWater
from oil_water_filtration.Enums import Components
from oil_water_filtration.Write import write, write_nevyaz


class OilWaterSS(OilWater):
    def show_results(self, time, layer, cell_container):
        x = np.arange(layer.x_0, layer.x_N, layer.h)
        x = np.append(x, layer.x_N)
        file = open('s(x)_count_' + str(time) + '_days_ss_method.txt', 'w')
        for i in range(cell_container.get_len()):
            state = cell_container.get_cells()[i].get_cell_state_n_plus()
            file.write(str(x[i]) + ' ' + str(state.get_s_water()) + '\n')
        file.close()
        file = open('p(x)_count_' + str(time) + '_days_ss_method.txt', 'w')
        for i in range(cell_container.get_len()):
            state = cell_container.get_cells()[i].get_cell_state_n_plus()
            file.write(str(x[i]) + ' ' + str(state.get_pressure_oil()) + '\n')
        file.close()

    def update_saturation(self, cell_container, flow_array):
        for cell in cell_container.get_cells()[1:-1]:
            state_n_plus = cell.get_cell_state_n_plus()
            state_n_plus.set_pressure_cap(state_n_plus.get_pressure_oil() - state_n_plus.get_pressure_water())
            s_water_graph = cell.layer.count_s_water_graph()
            s_w = cell.layer.count_s_water(s_water_graph, state_n_plus.get_pressure_cap())
            cell.get_cell_state_n_plus().set_s_water(s_w)
            cell.get_cell_state_n_plus().set_s_oil(1.0 - s_w)
            cells = cell_container.get_cells()
            cells[-1].get_cell_state_n_plus().set_s_water(cells[-2].get_cell_state_n_plus().get_s_water())
            cells[-1].get_cell_state_n_plus().set_s_oil(cells[-2].get_cell_state_n_plus().get_s_oil())

    @staticmethod
    def update_pressure(cell_container, delta_list_k):
        for i in range(1, cell_container.get_len() - 1):
            cell = cell_container.get_cells()[i]
            cell.get_cell_state_n_plus().set_pressure_water(
                cell.get_cell_state_n_plus().get_pressure_water() + delta_list_k[i - 1][0][0])
            cell.get_cell_state_n_plus().set_pressure_oil(
                cell.get_cell_state_n_plus().get_pressure_oil() + delta_list_k[i - 1][1][0])

    def count_max_elem(self, list_delta):
        list_delta_abs = []
        for elem in list_delta:
            list_delta_abs.append(math.fabs(elem))
        return max(list_delta_abs)

    #Условие на на оба давления
    def check_pressure_convergence(self, cell_container, delta_list_k):
        water_conv_list = []
        oil_conv_list = []
        for i in range(1, cell_container.get_len() - 1):
            cell = cell_container.get_cells()[i]
            water_conv = delta_list_k[i - 1][0][0] / cell.get_cell_state_n().get_pressure_water()
            oil_conv = delta_list_k[i - 1][1][0] / cell.get_cell_state_n().get_pressure_oil()
            water_conv_list.append(water_conv)
            oil_conv_list.append(oil_conv)

        if self.count_max_elem(water_conv_list) > 0.01 or self.count_max_elem(oil_conv_list) > 0.01:
            print("P_n+1 - P_n / P_n > 0.1")
            return True
        else:
            return False

    #Условие на капилярное давление
    def check_pressure_cap_convergence(self, cell_container, delta_list_k):
        p_cap_conv_list = []
        for i in range(1, cell_container.get_len() - 1):
            cell = cell_container.get_cells()[i]
            p_cap_conv = (delta_list_k[i - 1][1][0] - delta_list_k[i - 1][0][0]) / cell.get_cell_state_n().get_pressure_cap()
            p_cap_conv_list.append(p_cap_conv)

        if self.count_max_elem(p_cap_conv_list) > 4.0:
            print("P_n+1 - P_n / P_n > 0.1")
            return True
        else:
            return False

    @staticmethod
    def solve_slau(solver_slau):
        solver_slau.solve_thomas_matrix_method(0.0, 0.0, 0.0, 0.0)

    def generate_delta_k(self, layer):
        return [[0.0, 0.0]] + [[self.delta_0, self.delta_0] for i in range(layer.N - 2)] + [[0.0, 0.0]]

    @staticmethod
    def count_norm(list_delta):
        list_delta_abs = []
        for elem in list_delta:
            for el in elem:
                list_delta_abs.append(math.fabs(el))
        return max(list_delta_abs)

    @staticmethod
    def count_b(solverSlau, cell_flow):
        left_cell = cell_flow.get_left_cell()
        right_cell = cell_flow.get_right_cell()

        b = np.matrix(
            [[cell_flow.t_oil_water[Components.WATER.value], 0.0], [0.0, cell_flow.t_oil_water[Components.OIL.value]]])
        # от b
        t_p_11 = cell_flow.t_oil_water[Components.WATER.value] * right_cell.get_cell_state_n_plus().get_pressure_water()
        t_p_12 = cell_flow.t_oil_water[Components.OIL.value] * right_cell.get_cell_state_n_plus().get_pressure_oil()
        t_p = np.matrix([[t_p_11], [t_p_12]])

        # от а
        t_p_11_a = cell_flow.t_oil_water[
                       Components.WATER.value] * left_cell.get_cell_state_n_plus().get_pressure_water()
        t_p_12_a = cell_flow.t_oil_water[Components.OIL.value] * left_cell.get_cell_state_n_plus().get_pressure_oil()
        t_p_a = np.matrix([[t_p_11_a], [t_p_12_a]])

        solverSlau.set_matrix_coefficients_ss(cell_flow.get_left_cell().get_eq_index()[0],
                                              cell_flow.get_right_cell().get_eq_index()[0], b)
        if 0 <= left_cell.get_eq_index()[0] < solverSlau.e_count:

            solverSlau.add_nevyaz(left_cell.get_eq_index()[0], -t_p)
            solverSlau.add_nevyaz(left_cell.get_eq_index()[0], t_p_a)

    @staticmethod
    def count_c(solverSlau, cell_flow):
        right_cell = cell_flow.get_right_cell()
        left_cell = cell_flow.get_left_cell()
        c = np.matrix(
            [[cell_flow.t_oil_water[Components.WATER.value], 0], [0, cell_flow.t_oil_water[Components.OIL.value]]])

        # от с
        t_p_11 = cell_flow.t_oil_water[Components.WATER.value] * left_cell.get_cell_state_n_plus().get_pressure_water()
        t_p_12 = cell_flow.t_oil_water[Components.OIL.value] * left_cell.get_cell_state_n_plus().get_pressure_oil()
        t_p = np.matrix([[t_p_11], [t_p_12]])

        # от а
        t_p_11_a = cell_flow.t_oil_water[Components.WATER.value] * right_cell.get_cell_state_n_plus().get_pressure_water()
        t_p_12_a = cell_flow.t_oil_water[Components.OIL.value] * right_cell.get_cell_state_n_plus().get_pressure_oil()
        t_p_a = np.matrix([[t_p_11_a], [t_p_12_a]])

        solverSlau.set_matrix_coefficients_ss(cell_flow.get_right_cell().get_eq_index()[0],
                                              cell_flow.get_left_cell().get_eq_index()[0], c)
        if 0 <= right_cell.get_eq_index()[0] < solverSlau.e_count:
            solverSlau.add_nevyaz(right_cell.get_eq_index()[0], -t_p)
            solverSlau.add_nevyaz(right_cell.get_eq_index()[0], t_p_a)

    def count_d_left_diag(self, cell, s_component_n, ro_component_der, ro_component_n_plus):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        layer = cell.layer
        p_cap_graph = layer.count_pcap_graph()
        s_der = layer.count_s_water_graph_der(p_cap_graph, state_n.get_s_water())

        d = (1.0 / self.tau) * ((state_n.get_fi() * s_component_n) * ro_component_der
                                       - (state_n_plus.get_fi() * ro_component_n_plus) * s_der +
                                       0.5 * (ro_component_n_plus * s_component_n * layer.fi_0 * layer.c_r))
        return d

    def count_d_right_diag(self, cell, s_component_n, ro_component_n_plus):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        layer = cell.layer
        p_cap_graph = layer.count_pcap_graph()
        s_der = layer.count_s_water_graph_der(p_cap_graph, state_n.get_s_water())

        d = (1.0 / self.tau) * ((state_n_plus.get_fi() * ro_component_n_plus) * s_der
                                       + 0.5 * (ro_component_n_plus * s_component_n * layer.fi_0 * layer.c_r))
        return d

    def count_a_ss(self, solverSlau, cell):
        s_water_n = cell.get_cell_state_n().get_s_water()
        ro_water_der = cell.layer.ro_water_0 * cell.layer.c_f_water
        ro_water_n_plus = cell.get_cell_state_n_plus().get_ro_water()

        s_oil_n = cell.get_cell_state_n().get_s_oil()
        ro_oil_der = cell.layer.ro_oil_0 * cell.layer.c_f_oil
        ro_oil_n_plus = cell.get_cell_state_n_plus().get_ro_oil()

        d_11 = self.count_d_left_diag(cell, s_water_n, ro_water_der, ro_water_n_plus)
        d_12 = self.count_d_right_diag(cell, s_water_n, ro_water_n_plus)
        d_21 = self.count_d_right_diag(cell, s_oil_n, ro_oil_n_plus)
        d_22 = self.count_d_left_diag(cell, s_oil_n, ro_oil_der, ro_oil_n_plus)
        d = np.matrix([[d_11, d_12], [d_21, d_22]])
        d_p_11 = d_11 * (
            cell.get_cell_state_n_plus().get_pressure_water() - cell.get_cell_state_n().get_pressure_water()) \
                 + d_12 * (cell.get_cell_state_n_plus().get_pressure_oil() - cell.get_cell_state_n().get_pressure_oil())
        d_p_12 = d_21 * (
            cell.get_cell_state_n_plus().get_pressure_water() - cell.get_cell_state_n().get_pressure_water()) \
                 + d_22 * (cell.get_cell_state_n_plus().get_pressure_oil() - cell.get_cell_state_n().get_pressure_oil())
        d_p = np.matrix([[d_p_11], [d_p_12]])

        if 0 <= cell.get_eq_index()[0] < solverSlau.e_count:
            solverSlau.add_coefficient(cell.get_eq_index()[0], cell.get_eq_index()[0], -d)
            solverSlau.add_nevyaz(cell.get_eq_index()[0], d_p)

    def give_result_to_cells(self):
        pass

    def generate_matrix(self, flow_array, cell_container, solver_slau):
        for flow in flow_array:
            self.count_b(solver_slau, flow)
            self.count_c(solver_slau, flow)

        coeff_mat = write(solver_slau.coefficient_matrix)
        nevyaz_mat = write_nevyaz(solver_slau.nevyaz_vector)

        for cell in cell_container.get_cells():
            self.count_a_ss(solver_slau, cell)

        coeff_mat = write(solver_slau.coefficient_matrix)
        nevyaz_mat = write_nevyaz(solver_slau.nevyaz_vector)


    @staticmethod
    def reshape_matrices(solver_slau):
        solver_slau.coefficient_matrix = np.zeros((98, 98, 2, 2), float)
        solver_slau.nevyaz_vector = np.zeros((98, 2, 1))


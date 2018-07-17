import numpy as np
import math
from oil_water_filtration.OilWater import OilWater
from oil_water_filtration.Enums import Components


class OilWaterIMPES(OilWater):
    @staticmethod
    def update_pressure_cap(cell_container):
        for cell in cell_container.get_cells()[1:-1]:
            s_water = cell.get_cell_state_n_plus().get_s_water()
            #if s_water > 1.0:
                #s_water = 1.0
            p_cap_graph = cell.layer.count_pcap_graph()
            p_cap = cell.layer.count_pressure_cap(p_cap_graph, s_water)
            cell.get_cell_state_n_plus().set_pressure_cap(p_cap)


    def check_saturation_convergence(self, cell_container):
        if self.count_norm(
                [cell.get_cell_state_n_plus().get_s_oil() - cell.get_cell_state_n().get_s_oil() for cell in
                 cell_container.get_cells()]) > 0.1:
            print("S_n+1 - S_n > 0.1")
            return True
        else:
            return False

    def check_saturation_sum(self, cell_container):
        norm_list = [cell.get_cell_state_n_plus().get_s_water() + cell.get_cell_state_n_plus().get_s_oil() - 1.0 for cell in
                 cell_container.get_cells()]
        norm = self.count_norm(
                [cell.get_cell_state_n_plus().get_s_water() + cell.get_cell_state_n_plus().get_s_oil() - 1.0 for cell in
                 cell_container.get_cells()])
        if norm > 10 ** (-3):
            print("S_n_oil + S_n_water != 1     ","NORM    ", norm)
            return True
        else:
            return False

    def check_pressure_convergence(self, cell_container, delta_list_k):
        oil_conv_list = []
        for i in range(1, cell_container.get_len() - 1):
            cell = cell_container.get_cells()[i]
            oil_conv = delta_list_k[i] / cell.get_cell_state_n().get_pressure_oil()
            oil_conv_list.append(oil_conv)

        if self.count_norm(oil_conv_list) > 0.1:
            print("P_n+1 - P_n / P_n > 0.1")
            return True
        else:
            return False

    @staticmethod
    def solve_slau(solver_slau):
        solver_slau.solve_thomas_method(0.0, 0.0)

    def generate_delta_k(self, layer):
        return [0.0] + [self.delta_0 for i in range(layer.N - 2)] + [0.0]

    @staticmethod
    def count_norm(list_delta):
        list_delta_abs = [math.fabs(elem) for elem in list_delta[1:-1]]
        return max(list_delta_abs)

    @staticmethod
    def update_pressure(cell_container, delta_list_k):
        for i in range(1, cell_container.get_len()):
            cell = cell_container.get_cells()[i]
            state = cell.get_cell_state_n_plus()
            state.set_pressure_oil(state.get_pressure_oil() + delta_list_k[i])
            state.set_pressure_water(state.get_pressure_oil() - state.get_pressure_cap())




    @staticmethod
    #!!!!!!!!
    def count_a(solverSlau, cell):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        p_cap_graph = cell.layer.count_pcap_graph()
        p_cap_der = cell.layer.count_p_cap_graph_der(p_cap_graph, state_n.get_s_water())
        fi_der = cell.layer.fi_0 * cell.layer.c_r
        ro_der = cell.layer.ro_water_0 * cell.layer.c_f_water
        A = (state_n_plus.get_fi() * state_n_plus.get_ro_oil()) / (state_n_plus.get_fi() * state_n_plus.get_ro_water() - ro_der * state_n.get_s_water() * state_n.get_fi() * p_cap_der)
        solverSlau.set_matrix_coeffitients_a_d(cell.get_eq_index()[0], cell.get_eq_index()[0], state_n_plus.get_pressure_oil(), state_n.get_pressure_oil(), A * state_n_plus.get_c1_p() + state_n_plus.get_c2_p())


    @staticmethod
    #!!!!!!!!
    def count_b(solverSlau, cell_flow):
        left_cell = cell_flow.get_left_cell()
        right_cell = cell_flow.get_right_cell()
        state_n_plus = left_cell.get_cell_state_n_plus()
        state_n = left_cell.get_cell_state_n()
        p_cap_graph = left_cell.layer.count_pcap_graph()
        p_cap_der = left_cell.layer.count_p_cap_graph_der(p_cap_graph, state_n.get_s_water())
        fi_der = left_cell.layer.fi_0 * left_cell.layer.c_r
        ro_der = left_cell.layer.ro_water_0 * left_cell.layer.c_f_water
        A = (state_n_plus.get_fi() * state_n_plus.get_ro_oil()) / (state_n_plus.get_fi() * state_n_plus.get_ro_water() - ro_der * state_n.get_s_water() * state_n.get_fi() * p_cap_der)
        b = A * cell_flow.t_oil_water[Components.WATER.value] + cell_flow.t_oil_water[Components.OIL.value]
        r_ost = - A * cell_flow.t_oil_water[Components.WATER.value] * (right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        solverSlau.add_nevyaz(left_cell.get_eq_index()[0], -r_ost)
        solverSlau.set_matrix_coefficients(cell_flow.get_left_cell().get_eq_index()[0], cell_flow.get_right_cell().get_eq_index()[0], left_cell.get_cell_state_n_plus().get_pressure_oil(), right_cell.get_cell_state_n_plus().get_pressure_oil(), b)
        return b

    @staticmethod
    #!!!!!!!!!!!
    def count_c(solverSlau, cell_flow):
        left_cell = cell_flow.get_left_cell()
        right_cell = cell_flow.get_right_cell()
        state_n_plus = right_cell.get_cell_state_n_plus()
        state_n = right_cell.get_cell_state_n()
        p_cap_graph = right_cell.layer.count_pcap_graph()
        p_cap_der = right_cell.layer.count_p_cap_graph_der(p_cap_graph, state_n.get_s_water())
        fi_der = right_cell.layer.fi_0 * right_cell.layer.c_r
        ro_der = right_cell.layer.ro_water_0 * right_cell.layer.c_f_water
        A = (state_n_plus.get_fi() * state_n_plus.get_ro_oil()) / (state_n_plus.get_fi() * state_n_plus.get_ro_water() - ro_der * state_n.get_s_water() * state_n.get_fi() * p_cap_der)
        c = A * cell_flow.t_oil_water[Components.WATER.value] + cell_flow.t_oil_water[Components.OIL.value]
        r_ost = - A * cell_flow.t_oil_water[Components.WATER.value] * (
                    left_cell.get_cell_state_n_plus().get_pressure_cap() - right_cell.get_cell_state_n_plus().get_pressure_cap())
        solverSlau.add_nevyaz(right_cell.get_eq_index()[0], -r_ost)
        solverSlau.set_matrix_coefficients(cell_flow.get_right_cell().get_eq_index()[0],
                                           cell_flow.get_left_cell().get_eq_index()[0],
                                           right_cell.get_cell_state_n_plus().get_pressure_oil(),
                                           left_cell.get_cell_state_n_plus().get_pressure_oil(), c)
        return c

    def update_saturation(self, cell_container, flow_array):
        for i in range(1, cell_container.get_len() - 1):
            cells = cell_container.get_cells()
            cell_i = cells[i]
            cell_i_plus = cells[i + 1]
            cell_i_minus = cells[i - 1]

            cell_state_n_plus_i = cell_i.get_cell_state_n_plus()
            cell_state_n_plus_i_plus = cell_i_plus.get_cell_state_n_plus()
            cell_state_n_plus_i_minus = cell_i_minus.get_cell_state_n_plus()
            cell_state_n_i = cell_i.get_cell_state_n()

            #water
            ro_der = cell_i.layer.ro_water_0 * cell_i.layer.c_f_water
            fi_der = cell_i.layer.fi_0 * cell_i.layer.c_r

            #oil
            ro_der_oil = cell_i.layer.ro_oil_0 * cell_i.layer.c_f_oil

            p_cap_graph = cell_i.layer.count_pcap_graph()
            p_cap_der = cell_i.layer.count_p_cap_graph_der(p_cap_graph, cell_state_n_i.get_s_water())
            #p_cap_der = 1.0 / (cell_i.layer.count_s_water_graph_der(p_cap_graph, cell_state_n_i.get_s_water()))

            d_11 = (1.0 / self.tau) * cell_state_n_i.get_s_water() * (cell_state_n_i.get_fi() * ro_der + fi_der * cell_state_n_plus_i.get_ro_water())
            d_12 = (1.0 / self.tau) * (cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_water() - ro_der * cell_state_n_i.get_s_water() * cell_state_n_i.get_ro_water() * p_cap_der)
            d_21 = (1.0 / self.tau) * (1.0 - cell_state_n_plus_i.get_s_water()) * (cell_state_n_i.get_fi() * ro_der_oil + fi_der * cell_state_n_plus_i.get_ro_oil())
            d_22 = (1.0 / self.tau) * (cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_oil())

            coeff2 = flow_array[i].t_oil_water[Components.WATER.value] * (
                        cell_state_n_plus_i_plus.get_pressure_oil() - cell_state_n_plus_i.get_pressure_oil() + cell_state_n_plus_i.get_pressure_cap() - cell_state_n_plus_i_plus.get_pressure_cap()) + \
                     flow_array[i - 1].t_oil_water[Components.WATER.value] * (
                                 cell_state_n_plus_i_minus.get_pressure_oil() - cell_state_n_plus_i.get_pressure_oil() + cell_state_n_plus_i.get_pressure_cap() - cell_state_n_plus_i_minus.get_pressure_cap())  # !!!!!!!!!!!!!!!!!!!!!!add to t_water

            coeff_oil_check = flow_array[i].t_oil_water[Components.OIL.value] * (
                        cell_state_n_plus_i_plus.get_pressure_oil() - cell_state_n_plus_i.get_pressure_oil()) + \
                     flow_array[i - 1].t_oil_water[Components.OIL.value] * (
                                 cell_state_n_plus_i_minus.get_pressure_oil() - cell_state_n_plus_i.get_pressure_oil())

            coeff3_water = -(d_11 / d_12) * (cell_state_n_plus_i.get_pressure_oil() - cell_state_n_i.get_pressure_oil())

            coeff3_oil = (d_21 / d_22) * (cell_state_n_plus_i.get_pressure_oil() - cell_state_n_i.get_pressure_oil())

            #s_water_new = cell_state_n_i.get_s_water() + (
             #           self.tau / (cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_water())) * coeff2 + coeff3_water
            s_water_new = cell_state_n_i.get_s_water() + (1.0 / (d_12)) * coeff2 + coeff3_water

            s_oil_check = cell_state_n_i.get_s_oil() + (
                        self.tau / (cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_oil())) * coeff_oil_check + coeff3_oil

            #s_oil_check_list.append(s_oil_check)

            #cell_state_n_plus_i.set_s_oil(s_oil_check)
            #s_water_new = 1.0 - s_oil_check
            #if s_water_new < 10 ** (-6):
            #    s_water_new = 0.0
            #cell_state_n_plus_i.set_s_water(s_water_new)

            cell_state_n_plus_i.set_s_water(s_water_new)
            cell_state_n_plus_i.set_s_oil(1.0 - s_water_new)
            #cell_state_n_plus_i.set_s_oil(s_oil_check)
            #cell_state_n_plus_i.set_s_water(1.0 - s_oil_check)
            cells = cell_container.get_cells()
            cells[-1].get_cell_state_n_plus().set_s_water(cells[-2].get_cell_state_n_plus().get_s_water())
            cells[-1].get_cell_state_n_plus().set_s_oil(cells[-2].get_cell_state_n_plus().get_s_oil())


    def give_result_to_cells(self):
        pass

    def generate_matrix(self, flow_array, cell_container, solver_slau):
        for flow in flow_array:
            self.count_b(solver_slau, flow)
            self.count_c(solver_slau, flow)
        for cell in cell_container.get_cells():
            cell.get_cell_state_n_plus().set_c1_p(cell.layer.count_c1_p_new(self, cell))
            cell.get_cell_state_n_plus().set_c2_p(cell.layer.count_c2_p_new(self, cell))
            self.count_a(solver_slau, cell)

    def show_results(self, time, layer, cell_container):
        x = np.arange(layer.x_0, layer.x_N, layer.h)
        x = np.append(x, layer.x_N)
        file = open('s(x)_count_364_days.txt', 'w')
        file_pressure = open('pressure_water_impes_364.txt', 'w')
        for i in range(cell_container.get_len()):
            state = cell_container.get_cells()[i].get_cell_state_n_plus()
            file.write(str(x[i]) + ' ' + str(state.get_s_water()) + '\n')
            file_pressure.write(str(x[i]) + ' ' + str(state.get_pressure_water()) + '\n')
        file.close()
        file_pressure.close()


import numpy as np
import math
from oil_water_filtration.OilWater import OilWater
from oil_water_filtration.Enums import Components


class OilWaterIMPES(OilWater):
    def check_saturation_convergence(self, cell_container):
        if self.count_norm(
                [cell.get_cell_state_n_plus().get_s_water() - cell.get_cell_state_n().get_s_water() for cell in
                 cell_container.get_cells()]) > 0.1:
            print("S_n+1 - S_n > 0.1")
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
        for i in range(1, cell_container.get_len() - 1):
            cell = cell_container.get_cells()[i]
            state = cell.get_cell_state_n_plus()
            state.set_pressure_oil(state.get_pressure_oil() + delta_list_k[i])
            state.set_pressure_water(state.get_pressure_oil() - state.get_pressure_cap())



    @staticmethod
    def count_a(solverSlau, cell):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        A = state_n_plus.get_ro_oil() / state_n_plus.get_ro_water()
        solverSlau.set_matrix_coeffitients_a_d(cell.get_eq_index()[0], cell.get_eq_index()[0], state_n_plus.get_pressure_oil(), state_n.get_pressure_oil(), A * state_n_plus.get_c1_p() + state_n_plus.get_c2_p())


    @staticmethod
    def count_b(solverSlau, cell_flow):
        left_cell = cell_flow.get_left_cell()
        right_cell = cell_flow.get_right_cell()
        state_n_plus = left_cell.get_cell_state_n_plus()
        A = state_n_plus.get_ro_oil() / state_n_plus.get_ro_water()
        b = A * cell_flow.t_oil_water[Components.WATER.value] + cell_flow.t_oil_water[Components.OIL.value]
        r_ost = - A * cell_flow.t_oil_water[Components.WATER.value] * (right_cell.get_cell_state_n_plus().get_pressure_cap() - left_cell.get_cell_state_n_plus().get_pressure_cap())
        solverSlau.add_nevyaz(left_cell.get_eq_index()[0], -r_ost)
        solverSlau.set_matrix_coefficients(cell_flow.get_left_cell().get_eq_index()[0], cell_flow.get_right_cell().get_eq_index()[0], left_cell.get_cell_state_n_plus().get_pressure_oil(), right_cell.get_cell_state_n_plus().get_pressure_oil(), b)
        return b

    @staticmethod
    def count_c(solverSlau, cell_flow):
        left_cell = cell_flow.get_left_cell()
        right_cell = cell_flow.get_right_cell()
        state_n_plus = right_cell.get_cell_state_n_plus()
        A = state_n_plus.get_ro_oil() / state_n_plus.get_ro_water()
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

            coeff1 = cell_state_n_i.get_fi() * cell_state_n_i.get_ro_water() / (
                        cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_water())
            coeff2 = flow_array[i].t_oil_water[Components.WATER.value] * (
                        cell_state_n_plus_i_plus.get_pressure_oil() - cell_state_n_plus_i.get_pressure_oil() + cell_state_n_plus_i.get_pressure_cap() - cell_state_n_plus_i_plus.get_pressure_cap()) + \
                     flow_array[i - 1].t_oil_water[Components.WATER.value] * (
                                 cell_state_n_plus_i_minus.get_pressure_oil() - cell_state_n_plus_i.get_pressure_oil() + cell_state_n_plus_i.get_pressure_cap() - cell_state_n_plus_i_minus.get_pressure_cap())  # !!!!!!!!!!!!!!!!!!!!!!add to t_water
            s_water_new = cell_state_n_i.get_s_water() * coeff1 + (
                        self.tau / (cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_water())) * coeff2
            cell_state_n_plus_i.set_s_water(s_water_new)
            cell_state_n_plus_i.set_s_oil(1.0 - s_water_new)
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


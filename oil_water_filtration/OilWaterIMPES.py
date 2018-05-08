from oil_water_filtration.OilWater import OilWater
from oil_water_filtration.Enums import Components
import numpy as np

class OilWaterIMPES(OilWater):
    @staticmethod
    def update_pressure(cell_container, delta_list_k):
        for i in range(1, cell_container.get_len()):
            cell = cell_container.get_cells()[i]
            state = cell.get_cell_state_n_plus()
            state.set_pressure_oil(state.get_pressure_oil() + delta_list_k[i])
            state.set_pressure_water(state.get_pressure_oil() - state.get_pressure_cap())

    @staticmethod
    def count_a(solverSlau, cell):
        #A = ro_oil[index] / ro_water[index]
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

    @staticmethod
    def update_saturation(cell_container, flow_array, solver):
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
                        solver.tau / (cell_state_n_plus_i.get_fi() * cell_state_n_plus_i.get_ro_water())) * coeff2
            cell_state_n_plus_i.set_s_water(s_water_new)
            cell_state_n_plus_i.set_s_oil(1.0 - s_water_new)


    def give_result_to_cells(self):
        pass


    def generate_matrix(self, flow_array, cell_container, solver_slau, solver):
        for flow in flow_array:
            self.count_b(solver_slau, flow)
            self.count_c(solver_slau, flow)
        for cell in cell_container.get_cells():
            cell.get_cell_state_n_plus().set_c1_p(cell.layer.count_c1_p_new(solver, cell))
            cell.get_cell_state_n_plus().set_c2_p(cell.layer.count_c2_p_new(solver, cell))
            self.count_a(solver_slau, cell)

    @staticmethod
    def main_cycle(self):

        pass

    @staticmethod
    def show_results(oil_water_impes, time, layer, cell_container):
        print(time)
        time += oil_water_impes.tau
        oil_water_impes.tau = max(oil_water_impes.tau * 2.0, oil_water_impes.tau_default)
        if 125 * 86400.0 - time - oil_water_impes.tau <= 0:
            oil_water_impes.tau = 125 * 86400.0 - time  # - solver.tau

        if abs(time - 86400.0 * 125) < 1.e-16:
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
            for i in range(cell_container.get_len()):
                state = cell_container.get_cells()[i].get_cell_state_n_plus()
                file.write(str(x[i]) + ' ' + str(state.get_s_water()) + '\n')

            file.close()


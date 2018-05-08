from oil_water_filtration.OilWater import OilWater
import numpy as np
from oil_water_filtration.Enums import Components

class OilWaterSS(OilWater):
    @staticmethod
    def count_b(solverSlau, cell_flow):
        left_cell = cell_flow.get_left_cell()
        right_cell = cell_flow.get_right_cell()

        b = np.matrix([[cell_flow.t_oil_water[Components.WATER.value], 0], [0, cell_flow.t_oil_water[Components.OIL.value]]])
        #от b
        t_p_11 = cell_flow.t_oil_water[Components.WATER.value] * right_cell.get_cell_state_n_plus().get_pressure_water()
        t_p_12 = cell_flow.t_oil_water[Components.OIL.value] * right_cell.get_cell_state_n_plus().get_pressure_oil()
        t_p = np.matrix([[t_p_11], [t_p_12]])

        #от а
        t_p_11_a = cell_flow.t_oil_water[Components.WATER.value] * left_cell.get_cell_state_n_plus().get_pressure_water()
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
        left_cell = cell_flow.get_right_cell()
        c = np.matrix([[cell_flow.t_oil_water[Components.WATER.value], 0], [0, cell_flow.t_oil_water[Components.OIL.value]]])

        #от с
        t_p_11 = cell_flow.t_oil_water[Components.WATER.value] * left_cell.get_cell_state_n_plus().get_pressure_water()
        t_p_12 = cell_flow.t_oil_water[Components.OIL.value] * left_cell.get_cell_state_n_plus().get_pressure_oil()
        t_p = np.matrix([[t_p_11], [t_p_12]])

        #от а
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
        d = (layer.V_ij / self.tau) * ((state_n.get_fi() * s_component_n) * ro_component_der
                   - (state_n_plus.get_fi() * ro_component_n_plus) * (
                               (state_n_plus.get_s_water() - state_n.get_s_water()) / (
                                   state_n_plus.get_pressure_cap() - state_n.get_pressure_cap() + 1.0)) +#!!!!!!!!!!!!!!!!!
                   0.5 * (ro_component_n_plus * s_component_n * layer.fi_0 * layer.c_r))
        return d

    def count_d_right_diag(self, cell, s_component_n, ro_component_n_plus):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        layer = cell.layer
        d = (layer.V_ij / self.tau) * ((state_n_plus.get_fi() * ro_component_n_plus) * (
                    (state_n_plus.get_s_water() - state_n.get_s_water()) / (
                        state_n_plus.get_pressure_cap() - state_n.get_pressure_cap() + 1.0))#!!!!!!!!!!!!!
                   + 0.5 * (ro_component_n_plus * s_component_n * layer.fi_0 * layer.c_r))
        return d

    def count_a_ss(self, solverSlau, cell, delta):
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
        d_p_11 = d_11 * (cell.get_cell_state_n_plus().get_pressure_water() - cell.get_cell_state_n().get_pressure_water())\
                 + d_12 * (cell.get_cell_state_n_plus().get_pressure_oil() - cell.get_cell_state_n().get_pressure_oil())
        d_p_12 = d_21 * (cell.get_cell_state_n_plus().get_pressure_water() - cell.get_cell_state_n().get_pressure_water())\
                 + d_22 * (cell.get_cell_state_n_plus().get_pressure_oil() - cell.get_cell_state_n().get_pressure_oil())
        d_p = np.matrix([[d_p_11], [d_p_12]])

        if 0 <= cell.get_eq_index()[0] < solverSlau.e_count:
            solverSlau.add_coefficient(cell.get_eq_index()[0], cell.get_eq_index()[0], -d)
            solverSlau.add_nevyaz(cell.get_eq_index()[0], d_p)


    def give_result_to_cells(self):
        pass

    def generate_matrix(self):
        pass
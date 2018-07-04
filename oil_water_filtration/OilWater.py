from abc import abstractmethod
from oil_water_filtration.Enums import Components
import math

class OilWater:
    def __init__(self):
        self.delta_0 = 1000.0
        self.delta_max = 10 ** (-3)
        self.tau_default = 86400.0  # s- time step
        self.tau = self.tau_default
        # time_max = 365.25 * 4 * tau
        self.time_max = 728 * 1 * 86400

    @staticmethod
    @abstractmethod
    def count_norm(list_delta):
        pass

    @staticmethod
    def recount_properties(cell_container):
        for i in range(cell_container.get_len()):
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

    @staticmethod
    def count_flows(flows_array):
        for flow in flows_array:
            flow.count_cells_flows()

    @staticmethod
    @abstractmethod
    def solve_slau(solver_slau):
        pass

    @abstractmethod
    def generate_matrix(self, flow_array, cell_container, solverSlau):
        pass

    @abstractmethod
    def give_result_to_cells(self):
        pass

    @staticmethod
    @abstractmethod
    def count_a(solverSlau, cell):
        pass


    @staticmethod
    @abstractmethod
    def count_b(solverSlau, cell_flow):
        pass

    @staticmethod
    @abstractmethod
    def count_c(solverSlau, cell_flow):
        pass

    @staticmethod
    @abstractmethod
    def update_pressure(cell_container, delta_list_k):
        pass

    @abstractmethod
    def update_saturation(self, cell_container, flow_array):
        pass

    @staticmethod
    @abstractmethod
    def solve_slau():
        pass

    def main_cycle(self, layer, cell_container, solver_slau, flow_array):
        time = self.tau
        counter = 1
        while time < self.time_max:

            delta_list_k = self.generate_delta_k(layer)

            for cell in cell_container.get_cells():
                cell.get_cell_state_n().set_equals_to(cell.get_cell_state_n_plus())
            while self.count_norm(delta_list_k) > self.delta_max:
                solver_slau.set_zero()
                self.recount_properties(cell_container)
                self.count_flows(flow_array)
                self.generate_matrix(flow_array, cell_container, solver_slau)

                self.solve_slau(solver_slau)
                delta_list_k = solver_slau.get_result()
                solver_slau.clear_result()
                self.update_pressure(cell_container, delta_list_k)

            self.recount_properties(cell_container)
            self.count_flows(flow_array)
            self.update_saturation(cell_container, flow_array)

            if self.check_saturation_convergence(cell_container):
                self.tau = self.tau / 2.0
                for cell in cell_container.get_cells():
                    cell.get_cell_state_n_plus().set_equals_to(cell.get_cell_state_n())

            else:
                print(time)
                time += self.tau
                self.tau = max(self.tau * 2.0, self.tau_default)
                if 125 * 86400.0 - time - self.tau <= 0:
                    self.tau = 125 * 86400.0 - time

                if abs(time - 86400.0 * 125) < 1.e-16:
                    self.show_results(time, layer, cell_container)
                counter += 1


    @abstractmethod
    def show_results(self, time, layer, cell_container):
        pass

    @abstractmethod
    def generate_delta_k(self, layer):
        pass

    @abstractmethod
    def check_saturation_convergence(self, cell_container):
        pass


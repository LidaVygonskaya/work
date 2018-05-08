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
        self.time_max = 365.25 * 1 * 86400

    @staticmethod
    def count_norm(list_delta):
        list_delta_abs = [math.fabs(elem) for elem in list_delta[1:-1]]
        return max(list_delta_abs)

    @staticmethod
    def count_norm_ss(list_delta):
        list_delta_abs = [math.fabs(elem) for elem in list_delta]
        return max(list_delta_abs)


    @staticmethod
    def count_matrix_norm(list_delta):
        list_delta_abs = []
        for elem in list_delta:
            for el in elem:
                list_delta_abs.append(math.fabs(el))
        return max(list_delta_abs)

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

    def solve_slau(self):
        return 0

    @abstractmethod
    def generate_matrix(self, flow_array, cell_container, solverSlau, solver):
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

    @staticmethod
    @abstractmethod
    def solve_slau():
        pass

    @staticmethod
    @abstractmethod
    def main_cycle(self):
        pass

    @staticmethod
    @abstractmethod
    def show_results(self, time, layer,cell_container):
        pass


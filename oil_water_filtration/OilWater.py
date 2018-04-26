from abc import abstractmethod


class OilWater:
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
        t_oil_water = [[], []]
        for flow in flows_array:
            t_oil_water.append(flow.count_cells_flows()[0])
            t_oil_water.append(flow.count_cells_flows()[1])
        return t_oil_water

        return 0

    def solve_slau(self):
        return 0

    @abstractmethod
    def generate_matrix(self):
        pass

    @abstractmethod
    def give_result_to_cells(self):
        pass

    @abstractmethod
    def count_a(self):
        pass

    @abstractmethod
    def count_b(self):
        pass

    @abstractmethod
    def count_c(self):
        pass

    @abstractmethod
    def count_d(self):
        pass

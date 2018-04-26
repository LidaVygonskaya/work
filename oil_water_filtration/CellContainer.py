import numpy as np
from oil_water_filtration.Cell import Cell
from oil_water_filtration.CellState import CellState


class CellContainer:
    container = []

    def __init__(self, cell_count, layer):
        self.cell_count = cell_count
        for i in range(cell_count):
            if i == 0 or i == cell_count - 1:
                self.container.append(Cell(i, layer, -1))
            else:
                self.container.append(Cell(i, layer, i - 1))

    def get_cells(self):
        return self.container

    def get_len(self):
        return len(self.container)

    def initialize_cells(self):
        for cell in self.get_cells():
            pressure_cap_graph = cell.layer.count_pcap_graph()
            cell_index = self.get_cells().index(cell)
            for state in cell.get_cell_states():
                if cell_index == 0:
                    # Заполняем n-ый слой
                    state.set_s_water(cell.layer.s_water_left)
                    state.set_s_oil(cell.layer.s_oil_left)
                    state.set_pressure_oil(cell.layer.pressure_oil_left)
                    state.set_pressure_cap(cell.layer.count_pressure_cap(pressure_cap_graph, cell.layer.s_water_left))
                    state.set_pressure_water(cell.layer.pressure_water_left)

                elif cell_index == self.get_len() - 1:
                    state.set_s_water(cell.layer.s_water_right)
                    state.set_s_oil(cell.layer.s_oil_right)
                    state.set_pressure_oil(cell.layer.pressure_oil_right)
                    state.set_pressure_cap(cell.layer.count_pressure_cap(pressure_cap_graph, cell.layer.s_water_right))
                    state.set_pressure_water(cell.layer.pressure_water_right)

                else:
                    state.set_s_oil(cell.layer.s_oil_init)
                    state.set_s_water(cell.layer.s_water_init)
                    state.set_pressure_oil(cell.layer.pressure_oil_init)
                    state.set_pressure_cap(cell.layer.count_pressure_cap(pressure_cap_graph, cell.layer.s_water_init))
                    state.set_pressure_water(state.get_pressure_oil() - state.get_pressure_cap())

                for component in cell.layer.components:
                    component_index = cell.layer.components.index(component)
                    state.set_k_r(component.count_k_r(state.get_components_saturation()[component_index]), component_index)
                    state.set_ro(component.count_ro(state.get_components_pressure()[component_index]), component_index)

                state.set_fi(cell.layer.count_fi((state.get_pressure_oil() + state.get_pressure_water()) / 2.0))







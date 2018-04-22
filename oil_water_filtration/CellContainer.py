import numpy as np
from oil_water_filtration.Cell import Cell
from oil_water_filtration.CellState import CellState


class CellContainer:
    container = []

    def __init__(self, cell_count, layer):
        self.cell_count = cell_count
        for i in range(cell_count):
            self.container.append(Cell(i, layer))

    def get_cells(self):
        return self.container

    def get_len(self):
        return len(self.container)

    def initialize_cells(self):
        for cell in self.get_cells():
            cell_index = self.get_cells().index(cell)
            if cell_index == 0:
                # Заполняем n-ый слой
                cell.get_cell_state_n().set_pressure_oil(cell.layer.pressure_oil_left)
                cell.get_cell_state_n().set_pressure_water(cell.layer.pressure_water_left)
                cell.get_cell_state_n().set_s_oil(cell.layer.s_oil_left)
                cell.get_cell_state_n().set_s_water(cell.layer.s_water_left)

            elif cell_index == self.get_len() - 1:
                cell.get_cell_state_n().set_pressure_oil(cell.layer.pressure_oil_init)
                cell.get_cell_state_n().set_pressure_water(cell.layer.pressure_water_init)
                cell.get_cell_state_n().set_s_oil(cell.layer.s_oil_init)
                cell.get_cell_state_n().set_s_water(cell.layer.s_water_init)

            else:
                cell.get_cell_state_n().set_pressure_oil(cell.layer.pressure_oil_right)
                cell.get_cell_state_n().set_pressure_water(cell.layer.pressure_water_right)
                cell.get_cell_state_n().set_s_oil(cell.layer.s_oil_right)
                cell.get_cell_state_n().set_s_water(cell.layer.s_water_right)
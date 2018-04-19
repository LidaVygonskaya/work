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

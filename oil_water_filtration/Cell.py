from oil_water_filtration.CellState import CellState
from oil_water_filtration.Layer import Layer

class Cell:
    def __init__(self, x_coordinate, layer):
        self.layer = layer
        self.x_coordinate = x_coordinate
        self.width = self.layer.h
        self.cell_states = [CellState(), CellState()]#n, n + 1 layers

    #Cell constant parameters
    def get_k(self):
        return self.layer.k

    def get_c_f_water(self):
        return self.layer.c_f_water

    def get_c_f_oil(self):
        return self.layer.c_f_oil

    def get_mu_oil(self):
        return self.layer.mu_oil

    def get_mu_water(self):
        return self.layer.mu_water

    def get_mu_oil_water(self):
        return self.layer.mu_oil_water

    def get_cell_state_n(self):
        return self.cell_states[0]

    def get_cell_state_n_plus(self):
        return self.cell_states[1]














class Flow:
    def __init__(self, cell_1, cell_2):
        self.cell_1 = cell_1
        self.cell_2 = cell_2

    def get_max_pressure_cell(self, index):
        if self.cell_1.cell_state.pressure_oil_water[index] > self.cell_2.cell_state.pressure_oil_water[index]:
            return self.cell_1
        else:
            return self.cell_2

    def count_flow(self):
        t_oil_water = []
        for i in range(self.cell_1.layer.components_count):
            cell = self.get_max_pressure_cell(i)
            t_x = (cell.get_k() * cell.cell_states.get_k_r_oil_water[i] / cell.get_mu_oil_water[i]) * (1 / cell.layer.h) ** 2.0 * cell.get_ro_oil_water[i]
            t_oil_water.append(t_x)
        return t_oil_water




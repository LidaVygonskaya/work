from oil_water_filtration.Enums import Components


class Flow:
    def __init__(self, cell_1, cell_2):
        self.cell_1 = cell_1
        self.cell_2 = cell_2
        self.t_oil_water = [0.0, 0.0]
        self.t_oil_water_der = [0.0, 0.0]# Производная потока по капиллярному давлению

    def get_max_pressure_cell(self, index):
        cell_1_pressure = self.cell_1.get_cell_state_n_plus().get_components_pressure()
        cell_2_pressure = self.cell_2.get_cell_state_n_plus().get_components_pressure()
        if cell_1_pressure[index] > cell_2_pressure[index]:
            return self.cell_1
        else:
            return self.cell_2

    def count_cells_flows(self):
        for i in range(self.cell_1.layer.components_count):
            cell = self.get_max_pressure_cell(i)
            cell_state_n_plus = cell.get_cell_state_n_plus()
            cell_state_n = cell.get_cell_state_n()
            p_cap_graph = cell.layer.count_pcap_graph()
            s_w = cell_state_n_plus.get_s_water()
            s_der = cell.layer.count_s_water_graph_der(p_cap_graph, s_w)
            if i == Components.OIL.value:
                s_der = (-1.0) * s_der
            t_x = (cell.get_k() * cell_state_n_plus.get_components_k_r()[i] / cell.get_mu_oil_water()[i]) * (1 / cell.layer.h) ** 2.0 * cell_state_n_plus.get_components_ro()[i]
            t_x_der = 2.0 * s_der * (cell.get_k() / cell.get_mu_oil_water()[i]) * ((1 / cell.layer.h) ** 2.0) * cell_state_n_plus.get_components_ro()[i] * cell_state_n_plus.get_components_saturation()[i]
            self.t_oil_water[i] = t_x
            self.t_oil_water_der[i] = t_x_der

    def get_right_cell(self):
        return self.cell_2

    def get_left_cell(self):
        return self.cell_1

    def count_flows(self, cell_container):
        cells = cell_container.get_cells()
        for i in range(cell_container.get_len() - 1):
            flow = Flow(cells[i], cells[i + 1])
            self.count_cells_flows(flow)



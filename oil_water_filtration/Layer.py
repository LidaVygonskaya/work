from oil_water_filtration.Oil import Oil
from oil_water_filtration.Water import Water
from oil_water_filtration.Enums import Components

class Layer:
    def __init__(self):
        self.atm = 101325.0
        self.ro_oil_0 = 900.0
        self.ro_water_0 = 1000.0
        self.P_01 = 80 * self.atm  # new
        self.P_02 = 80 * self.atm
        self.c_f_oil = (10.0 ** (-4)) / self.atm
        self.c_f_water = (10.0 ** (-4)) / self.atm
        self.c_r = (10.0 ** (-5)) / self.atm
        self.fi_0 = 0.2

        self.mu_oil = 10.0 * (10.0 ** (-3))
        self.mu_water = 10.0 ** (-3)
        self.mu_oil_water = [self.mu_oil, self.mu_water]

        self.k = (9.868233 * (10 ** (-13))) * 10 ** (0)
        self.s_water_init = 0.0 #10 ** (-4)
        self.s_oil_init = 1.0 - self.s_water_init
        self.pressure_cap_init = [] #list
        self.pressure_oil_init = 80.0 * self.atm#!!!!!!!!!!!!!!
        self.pressure_water_init = 80 * self.atm
        self.components = [Oil(), Water()]

        self.components_count = len(self.components)
        #for m in range(len(self.pressure_cap_init[1: -1])):
        #    self.pressure_water_init.append(self.pressure_oil_init - self.pressure_cap_init[1: -1][m])



        #left side
        self.pressure_water_left = 130.0 * self.atm#!!!!!!!!!!!!!!!
        self.pressure_oil_left = 130.0 * self.atm
        self.s_water_left = 1.0
        self.s_oil_left = 0.0
        #self.pressure_oil_left = self.pressure_water_left + self.pressure_cap_init[0]

        #right side
        self.pressure_water_right = 80.0 * self.atm
        self.pressure_oil_right = 80.0 * self.atm
        self.s_oil_right = self.s_oil_init
        self.s_water_right = self.s_water_init

        #space
        self.x_0 = 0.0
        self.x_N = 100.0  # meters
        self.N = 50
        self.h = (self.x_N - self.x_0) / (self.N - 1)
        self.z = 10.0
        self.V_ij = self.h ** 2.0 * self.z


    def get_water_component(self):
        return self.components[Components.WATER.value]

    def get_oil_component(self):
        return self.components[Components.OIL.value]


    @staticmethod
    def count_k_r(_s_water):
        k_r_water = _s_water ** 2.0
        k_r_oil = (1.0 - _s_water) ** 2.0
        return [k_r_water, k_r_oil]

    def count_ro_water_oil(self,_pressure_water, _pressure_oil):
        ro_water = self.ro_water_0 * (1.0 + self.c_f_water * (_pressure_water - self.P_02))
        ro_oil = self.ro_oil_0 * (1.0 + self.c_f_oil * (_pressure_oil - self.P_02))
        return [ro_water, ro_oil]

    def count_ro_water(self, _pressure_water):
        return self.ro_water_0 * (1.0 + self.c_f_water * (_pressure_water - self.P_02))

    def count_ro_oil(self, _pressure_oil):
        return self.ro_oil_0 * (1.0 + self.c_f_oil * (_pressure_oil - self.P_02))

    def count_fi(self, pressure_oil):
        return self.fi_0 * (1.0 + self.c_r * (pressure_oil - self.P_01))

    def count_pcap_graph(self):
        pressure_cap_graph = {}  # from graph {s_w: pressure_cap}
        file = open('Pcap(Sw).txt', 'r')
        for line in file.readlines():
            line = line.rstrip()
            s = line.split('\t')
            pressure_cap_graph.update({float(s[0]): float(s[1]) * self.atm})
        return pressure_cap_graph

    def count_s_water_graph(self):
        s_water_graph = {}
        file = open('Pcap(Sw).txt', 'r')
        for line in file.readlines():
            line = line.rstrip()
            s = line.split('\t')
            s_water_graph.update({float(s[1]) * self.atm: float(s[0])})
        return s_water_graph

    #Считаем производную для SS метода
    @staticmethod
    def count_s_water_graph_der(p_cap_graph, s_water):
        s_w_graph = list(p_cap_graph.keys())
        for i in range(len(s_w_graph)):
            if s_water <= s_w_graph[i]:
                s_der = (s_w_graph[i] - s_w_graph[i - 1]) / (p_cap_graph.get(s_w_graph[i]) - p_cap_graph.get(s_w_graph[i - 1]))
                break

        return s_der

    @staticmethod
    def count_p_cap_graph_der(p_cap_graph, s_water):
        s_w_graph = list(p_cap_graph.keys())
        for i in range(len(s_w_graph)):
            if s_water <= s_w_graph[i]:
                p_cap_der = (p_cap_graph.get(s_w_graph[i]) - p_cap_graph.get(s_w_graph[i - 1])) / (s_w_graph[i] - s_w_graph[i - 1])
                break
        return p_cap_der


    @staticmethod
    def count_s_water_graph_der_from_p(s_water_graph, p_cap):
        p_cap_graph = list(reversed(list(s_water_graph.keys())))
        for i in range(len(p_cap_graph)):
            if p_cap <= p_cap_graph[i]:
                s_i = s_water_graph.get(p_cap_graph[i])
                s_i_minus = s_water_graph.get(p_cap_graph[i - 1])
                p_cap_i = p_cap_graph[i]
                p_cap_minus = p_cap_graph[i - 1]
                s_der = (s_water_graph.get(p_cap_graph[i]) - s_water_graph.get(p_cap_graph[i - 1])) / (p_cap_graph[i] - p_cap_graph[i - 1])
                break
        return s_der

    @staticmethod
    def count_pressure_cap(p_cap_graph, s_water):
        s_w_graph = list(p_cap_graph.keys())
        for i in range(len(s_w_graph)):
            if s_water <= s_w_graph[i]:
                p_cap = p_cap_graph.get(s_w_graph[i - 1]) + (p_cap_graph.get(s_w_graph[i]) - p_cap_graph.get(s_w_graph[i - 1])) / (s_w_graph[i] - s_w_graph[i - 1]) * (s_water - s_w_graph[i - 1])
                break
        return p_cap

    def count_s_water(self, s_water_graph, p_cap):
        p_cap_graph = list(s_water_graph.keys())
        p_cap_graph.reverse()
        for i in range(len(p_cap_graph)):
            if p_cap <= p_cap_graph[i]:
                s_wat = s_water_graph.get(p_cap_graph[i - 1]) + (s_water_graph.get(p_cap_graph[i]) - s_water_graph.get(p_cap_graph[i - 1])) / (p_cap_graph[i] - p_cap_graph[i - 1]) * (p_cap - p_cap_graph[i - 1])
                break
            s_wat = 0.0
        return s_wat

    # @staticmethod
    # def count_pcap(p_cap_graph, s_water_list):
    #     p_cap_list = []
    #     s_w_graph = list(p_cap_graph.keys())
    #     for i in range(len(s_water_list)):
    #         for j in range(len(s_w_graph)):
    #             if s_water_list[i] <= s_w_graph[j]:
    #                 p_cap = p_cap_graph.get(s_w_graph[j - 1]) + (p_cap_graph.get(s_w_graph[j]) - p_cap_graph.get(s_w_graph[j - 1])) / (s_w_graph[j] - s_w_graph[j - 1]) * (s_water_list[i] - s_w_graph[j - 1])
    #                 # p_cap_list.append(p_cap)
    #                 p_cap_list.append(0.0)
    #                 break
    #     return p_cap_list

    @staticmethod
    def count_c1_p_new(solver, cell):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        cell_layer = cell.layer
        return ((state_n.get_fi() * state_n.get_s_water() * cell_layer.ro_water_0 * cell_layer.c_f_water) + (state_n.get_s_water() * state_n_plus.get_ro_water() * cell_layer.c_r * cell_layer.fi_0)) / solver.tau

    @staticmethod
    def count_c2_p_new(solver, cell):
        state_n = cell.get_cell_state_n()
        state_n_plus = cell.get_cell_state_n_plus()
        cell_layer = cell.layer
        return ((1.0 - state_n.get_s_water()) * (state_n.get_fi() * cell_layer.ro_oil_0 * cell_layer.c_f_oil + cell_layer.c_r * cell_layer.fi_0 * state_n_plus.get_ro_oil())) / solver.tau
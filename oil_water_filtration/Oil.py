class Oil:
    def __init__(self):
        self.ro_oil_0 = 900.0
        self.c_f_oil = (10.0 ** (-4)) / 101325.0
        self.P_02 = 80 * 101325.0
        self.pressure_oil_init = 80 * 101325.0
        self.pressure_oil_left = 130 * 101325.0

    def count_ro(self, pressure_oil):
        return self.ro_oil_0 * (1.0 + self.c_f_oil * (pressure_oil - self.P_02))

    @staticmethod
    def count_k_r(s_oil):
        #s_oil = max(0.0, s_oil)
        #_oil = min(s_oil, 1.0)
        return s_oil ** 2.0




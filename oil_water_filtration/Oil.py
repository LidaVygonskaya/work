class Oil:
    def __init__(self):
        self.ro_oil_0 = self.ro_oil_0 = 900.0
        self.c_f_oil = 0.0
        self.P_02 = 80 * 101325.0
        self.pressure_oil_init = 80 * 101325.0
        self.pressure_oil_left = 130 * 101325.0

    def count_ro(self, pressure_oil):
        return self.ro_oil_0 * (1.0 + self.c_f_oil * (pressure_oil - self.P_02))
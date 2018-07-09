class Water:
    def __init__(self):
        self.ro_water_0 = 1000.0
        self.c_f_water = (10.0 ** (-4)) / 101325.0
        self.P_01 = 80 * 101325.0

    def count_ro(self, pressure_water):
        return self.ro_water_0 * (1.0 + self.c_f_water * (pressure_water - self.P_01))

    @staticmethod
    def count_k_r(s_water):
        s_water = max(0.0, s_water)
        s_water = min(s_water, 1.0)
        return s_water ** 2.0
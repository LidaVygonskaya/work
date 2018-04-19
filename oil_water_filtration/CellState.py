class CellState:
    def __init__(self):
        self.pressure_oil = 0.0
        self.pressure_water = 0.0
        self.pressure_oil_water = [self.pressure_oil, self.pressure_water]
        self.s_water = 0.0
        self.s_oil = 0.0
        self.s_oil_water = [self.s_oil, self.s_water]
        self.ro_water = 0.0
        self.ro_oil = 0.0
        self.ro_oil_water = [self.ro_oil, self.ro_water]
        self.k_r_water = 0.0
        self.k_r_oil = 0.0
        self.k_r_oil_water = [self.k_r_oil, self.k_r_water]
        self.fi = 0.0

    def get_pressure_oil(self):
        return self.pressure_oil

    def get_pressure_water(self):
        return self.pressure_water

    def get_k_r_oil_water(self):
        return self.k_r_oil_water

    def get_s_water(self):
        return self.s_water

    def get_s_oil(self):
        return self.s_oil

    def get_ro_water(self):
        return self.ro_water

    def get_ro_oil(self):
        return self.ro_oil

    def get_ro_oil_water(self):
        return self.ro_oil_water

    def get_k_r_oil(self):
        return self.k_r_oil

    def get_k_r_water(self):
        return self.k_r_water

    def set_pressure_oil(self, pressure_oil):
        self.pressure_oil = pressure_oil

    def set_pressure_water(self, pressure_water):
        self.pressure_water = pressure_water

    def set_s_water(self, s_water):
        self.s_water = s_water

    def set_s_oil(self, s_oil):
        self.s_oil = s_oil

    def set_ro_water(self, ro_water):
        self.ro_water = ro_water

    def set_ro_oil(self, ro_oil):
        self.ro_oil = ro_oil

    def set_fi(self, fi):
        self.fi = fi




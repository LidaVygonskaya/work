class CellState:
    def __init__(self):
        self.pressure_oil_water = [0.0, 0.0]
        self.s_oil_water = [0.0, 0.0]
        self.ro_oil_water = [0.0, 0.0]
        self.k_r_oil_water = [0.0, 0.0]
        self.fi = 0.0

    #Getters
    def get_pressure_oil(self):
        return self.pressure_oil_water[0]

    def get_pressure_water(self):
        return self.pressure_oil_water[1]

    def get_components_pressure(self):
        return self.pressure_oil_water

    def get_components_saturation(self):
        return self.s_oil_water

    def get_s_water(self):
        return self.s_oil_water[1]

    def get_s_oil(self):
        return self.s_oil_water[0]

    def get_ro_water(self):
        return self.ro_oil_water[1]

    def get_ro_oil(self):
        return self.ro_oil_water[0]

    def get_components_ro(self):
        return self.ro_oil_water

    def get_components_k_r(self):
        return self.k_r_oil_water

    def get_k_r_oil(self):
        return self.k_r_oil_water[0]

    def get_k_r_water(self):
        return self.k_r_oil_water[1]

    #Setters
    def set_pressure_oil(self, pressure_oil):
        self.pressure_oil_water[0] = pressure_oil

    def set_pressure_water(self, pressure_water):
        self.pressure_oil_water[1] = pressure_water

    def set_s_water(self, s_water):
        self.s_oil_water[1] = s_water

    def set_s_oil(self, s_oil):
        self.s_oil_water[0] = s_oil

    def set_ro_water(self, ro_water):
        self.ro_oil_water[1] = ro_water

    def set_ro_oil(self, ro_oil):
        self.ro_oil_water[0] = ro_oil

    def set_k_r(self, k_r, index):
        self.k_r_oil_water[index] = k_r

    def set_ro(self, ro, index):
        self.ro_oil_water[index] = ro

    def set_fi(self, fi):
        self.fi = fi


    #Set n to n plus 1
    def set_equals_to(self, cell_state_new):
        self.pressure_oil_water[0] = cell_state_new.pressure_oil_water[0]
        self.pressure_oil_water[1] = cell_state_new.pressure_oil_water[1]
        self.s_oil_water[0] = cell_state_new.s_oil_water[0]
        self.s_oil_water[1] = cell_state_new.s_oil_water[1]
        self.ro_oil_water[0] = cell_state_new.ro_oil_water[0]
        self.ro_oil_water[1] = cell_state_new.ro_oil_water[1]
        self.k_r_oil_water[0] = cell_state_new.k_r_oil_water[0]
        self.k_r_oil_water[1] = cell_state_new.k_r_oil_water[1]
        self.fi = cell_state_new.fi






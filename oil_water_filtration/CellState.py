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

    #Getters
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

    #Setters
    def set_pressure_oil(self, pressure_oil):
        self.pressure_oil = pressure_oil
        self.pressure_oil_water[0] = pressure_oil

    def set_pressure_water(self, pressure_water):
        self.pressure_water = pressure_water
        self.pressure_oil_water[1] = pressure_water

    def set_s_water(self, s_water):
        self.s_water = s_water
        self.s_oil_water[1] = s_water

    def set_s_oil(self, s_oil):
        self.s_oil = s_oil
        self.s_oil_water[0] = s_oil

    def set_ro_water(self, ro_water):
        self.ro_water = ro_water
        self.ro_oil_water[1] = ro_water

    def set_ro_oil(self, ro_oil):
        self.ro_oil = ro_oil
        self.ro_oil_water[0] = ro_oil

    def set_fi(self, fi):
        self.fi = fi


    #Set n to n plus 1
    def set_equals_to(self, cell_state_new):
        self.pressure_oil = cell_state_new.presssure_oil
        self.pressure_oil_water[0] = cell_state_new.presssure_oil_water[0]

        self.pressure_water = cell_state_new.pressure_water
        self.pressure_oil_water[1] = cell_state_new.pressure_oil_water[1]

        self.s_oil = cell_state_new.s_oil
        self.s_oil_water[0] = cell_state_new.s_oil_water[0]

        self.s_water = cell_state_new.s_water
        self.s_oil_water[1] = cell_state_new.s_oil_water[1]

        self.ro_oil = cell_state_new.ro_oil
        self.ro_oil_water[0] = cell_state_new.ro_oil_water[0]

        self.ro_water = cell_state_new.ro_water
        self.ro_oil_water[1] = cell_state_new.ro_oil_water[1]





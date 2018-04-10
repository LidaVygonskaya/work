from oil_water_filtration.Layer import Layer
import math
import numpy as np
import matplotlib.pyplot as plt


class AnalHelper:
    def __init__(self):
        self.layer = Layer()
        self.eta_0 = (self.layer.mu_water / self.layer.mu_oil) * (self.layer.ro_oil_0 / self.layer.ro_water_0)
        self.L = self.layer.x_N
        self.delta_p = 50.0 * self.layer.atm
        self.time = 125.0 * 86400.0

    def count_s_c(self):
        #return math.sqrt(self.eta_0 / (self.eta_0 + self.layer.ro_oil_0))
        return math.sqrt(self.eta_0 / (1.0 + self.eta_0))
        #return math.sqrt(self.layer.mu_water / (self.layer.mu_water + self.layer.mu_oil))

    def count_f_der(self, s):
        return (2.0 * self.eta_0 * s * (1.0 - s)) / ((s ** 2.0 + self.eta_0 * (1.0 - s) ** 2.0) ** 2.0)

    def integral_function(self, s):
        return (2.0 * self.count_f_der(s) * (s + self.eta_0 * (s - 1.0))) /((s ** 2.0 + self.eta_0 * (1 - s) ** 2.0) ** 2.0)

    def count_v_m_omega(self, psi_c, s_c, time):
        return (-self.L + math.sqrt(self.L ** 2.0 + (2.0 / self.layer.fi_0) * (self.eta_0 * psi_c - self.count_f_der(s_c)) * self.layer.k * self.delta_p * time / self.layer.mu_oil)) / (self.eta_0 * psi_c - self.count_f_der(s_c))

    def count_x_c(self, s_c, v_m_omega):
        return v_m_omega * self.count_f_der(s_c)

    def count_x(self, V_m_omega, s):
        return self.count_f_der(s) * V_m_omega

    def count_integral(self, s_c):
        step = (s_c - 1.0) / self.layer.N
        s_grid = np.arange(1.0, s_c, step)
        s_grid = np.append(s_grid, [s_c])
        summ = 0.0
        integral = 0.0

        for i in range(self.layer.N):
            summ += 0.5 * (self.integral_function(s_grid[i]) + self.integral_function(s_grid[i + 1])) * (s_grid[i + 1] - s_grid[i])
        integral += summ
        return integral


analHelper = AnalHelper()
s_c = analHelper.count_s_c()
psi_c = analHelper.count_integral(s_c) + analHelper.count_f_der(s_c) / (s_c ** 2.0 + analHelper.eta_0 * (1.0 - s_c) ** 2.0)
V_m_omega = analHelper.count_v_m_omega(psi_c, s_c, analHelper.time)
step = (s_c - 1.0) / 1000.0
s_grid = np.arange(1.0, s_c, step)
s_grid = np.append(s_grid, [s_c])
x_list = [analHelper.count_x(V_m_omega, s_grid[i]) for i in range(len(s_grid))]
file = open('s(x)_anal_125_days.txt', 'w')
for i in range(len(x_list)):
    file.write(str(x_list[i]) + '   ' + str(s_grid[i]) + '\n')
file.close()
plt.plot(x_list, s_grid)
plt.show()
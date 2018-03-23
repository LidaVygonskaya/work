from matplotlib import cm
import two_dim
import numpy as np
from matplotlib import pyplot as plt


class ImplicitSolver:
    def __init__(self):
        self.tau = 720.0  # s- time step
        # time_max = 365.25 * 4 * tau
        self.time_max = 365.25 * 1 * 86400.0
        self.t_out_fields = 10 * self.tau  # frecuency of writing

    #@staticmethod
    def count_pressure(self, layer, a, b, c, f, g, pressure, beta_matrix,wi_matrix, well):
        pressure_n = np.ones((N_y + 1, N_x + 1)) * layer.P_init
        for i in range(layer.N_y - 1):
            for j in range(layer.N_x - 1):
                pressure_n[i + 1][j + 1] = (self.tau / (layer.V_ij * beta_matrix[i][j])) * (c[i][j] * pressure[i + 1][j] + g[i][j] * pressure[i][j + 1] + a[i][j] * pressure[i + 1][j + 1] + b[i][j] * pressure[i + 1][j + 2] + f[i][j] * pressure[i + 2][j + 1] + q_matrix[i][j])
                #pressure_n[i + 1][j + 1] = (1.0 / ((layer.V_ij * beta_matrix[i][j] / self.tau) - wi_matrix[i][j])) * (c[i][j] * pressure[i + 1][j] + g[i][j] * pressure[i][j + 1] + a[i][j] * pressure[i + 1][j + 1] + b[i][j] * pressure[i + 1][j + 2] + f[i][j] * pressure[i + 2][j + 1] - wi_matrix[i][j] * well.pressure_w)
        return pressure_n

    @staticmethod
    def count_q(q_matrix, well, wi_matrix, p_matrix):
        for i in range(layer.N_y - 1):
            for j in range(layer.N_x - 1):
                q_matrix[i][j] = wi_matrix[i][j] * (p_matrix[i][j] - well.pressure_w)
        return q_matrix


    def count_tau(self, c, b, g, f, beta_matrix, layer):
        matrix = np.ones((layer.N_y - 1, layer.N_x - 1))
        for i in range(layer.N_y - 1):
            for j in range(layer.N_x - 1):
                matrix[i, j] = (1 / 4.0) * (layer.V_ij * beta_matrix[i, j]) * (1 / (c[i, j] + b[i, j] + g[i, j] + f[i, j]))
        matrix = np.absolute(matrix)
        tau = np.amax(matrix)
        return tau

layer = two_dim.Layer()
implicit_solver = ImplicitSolver()
solver = two_dim.Solver()
solver.tau = implicit_solver.tau
well = two_dim.Well()

N_x = layer.N_x
N_y = layer.N_y
p_matrix = np.ones((N_y + 1, N_x + 1)) * layer.P_init #размер 101 на 101
q_matrix = np.ones((N_y - 1, N_x - 1))
time = solver.tau
p_list = []
q_list = []
time_list = []
counter = 1
while time < solver.time_max:
    beta = layer.count_beta(p_matrix)
    wi = well.count_wi(p_matrix, layer)
    q_matrix = implicit_solver.count_q(q_matrix, well, wi, p_matrix)
    c, b = solver.count_tx(p_matrix, layer)
    g, f = solver.count_ty(p_matrix, layer)
    a = -(b + c + g + f) + layer.V_ij * beta / solver.tau
    #t = implicit_solver.count_tau(c, b, g, f, beta, layer)
    #print(t)
    p_matrix_new = implicit_solver.count_pressure(layer, a, b, c, f, g, p_matrix, beta, wi, well)
    p_matrix = p_matrix_new.copy()
    if counter == 1200:
        p_list = p_matrix[45, :]
        X = np.arange(layer.x_0 + layer.h / 2, layer.x_N, layer.h)
        file_p = open('pressure_10_days_720.txt', 'w')
        for i in range(len(p_list[1:-1])):
            file_p.write(str(X[i]) + '  ' + str(p_list[1:-1][i]) + '\n')
        file_p.close()



    #fig = plt.figure()
    #ax = fig.gca(projection='3d')

    # Make data.
    # X = np.arange(layer.x_0 + layer.h / 2, layer.x_N, layer.h)
    # Y = np.arange(layer.y_0 + layer.h / 2, layer.y_N, layer.h)
    # X, Y = np.meshgrid(X, Y)
    #
    # # Plot the surface.
    # surf = ax.plot_surface(X, Y, p_matrix[1:-1, 1:-1], cmap=cm.coolwarm,
    #                        linewidth=0, antialiased=False)
    # #plt.savefig(name)
    # #plt.clf()
    # plt.show()
    time_list.append(time)
    q_list.append(q_matrix[45, 45])
    file = open('q_time_720_step.txt', 'w')
    for i in range(len(q_list)):
        file.write(str(time_list[i]) + '  ' + str(q_list[i]) + '\n')
    file.close()
    print(str(time) + ' ' + str(q_matrix[45, 45]))
    time += solver.tau
    counter += 1


    # #d = wi * well.pressure_w - p_n * (layer.V_ij * beta / solver.tau)


# X = np.arange(layer.x_0 + layer.h / 2, layer.x_N, layer.h)
# file_p = open('pressure_10_days_720.txt', 'w')
# for i in range(len(p_list)):
#     file.write(str(X[i]) + '  ' + str(p_list[1:-1][i]) + '\n')
# file.close()
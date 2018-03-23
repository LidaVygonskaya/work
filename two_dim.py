import math
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from matplotlib import pyplot as plt
import copy
# coding: utf8




class Layer:
    #Класс пласта
    def __init__(self):
        self.atm = 101325.0
        self.fi = 0.2  # poristost
        self.k = (9.868233 * (10 ** (-13))) * 10 ** (-3)
        self.k_x = self.k
        self.k_y = self.k
        self.P_init = 80.0 * self.atm  # Pa
        self.P_0 = 130 * self.atm
        self.P_Nx = 30 * self.atm
        self.ro_0 = 1000  # kg/m3 new
        self.P_01 = 80 * self.atm  # new
        self.P_02 = 80 * self.atm  # new
        self.fi_0 = 0.2
        #fluid
        self.c_r = (10 ** (-5)) / self.atm  # new
        self.c_f = (10 ** (-4)) / self.atm  # szhimaemost'
        self.mu = 10 ** (-3)  # Pa * s
        #space
        self.x_0 = 0.0
        self.x_N = 500.0  # meters
        self.y_0 = 0.0
        self.y_N = 500.0
        self.N_x = 100
        self.N_y = 100
        self.N = 100.0
        self.h = (self.x_N - self.x_0) / (self.N_x - 1)
        self.z = 10.0
        self.V_ij = self.h ** 2.0 * self.z

    def count_ro(self, pressure_ij):
        return self.ro_0 * (1.0 + self.c_f * (pressure_ij - self.P_02))

    def count_fi_p(self, pressure_ij):
        return self.fi_0 * (1.0 + self.c_r * (pressure_ij - self.P_01))

    def count_beta(self, pressure_matrix):
        beta = np.zeros((self.N_y - 1, self.N_x - 1))
        for i in range(self.N_y - 1):
            for j in range(self.N_x - 1):
                beta[i][j] = self.fi_0 * self.c_r * self.count_ro(pressure_matrix[i][j]) + self.c_f * self.ro_0 * self.count_fi_p(pressure_matrix[i][j])
        return beta


class Well:
    #Класс скважины
    def __init__(self):
        self.atm = 101325.0
        self.r_w = 0.108  # meters
        self.pressure_w = 79.0 * self.atm
        self.x_well = 45
        self.y_well = 45

    @staticmethod
    def count_re(layer):
        return layer.h * math.exp(-math.pi / 2.0)

    def count_wi(self, pressure_matrix, layer):
        wi = np.zeros((layer.N_x - 1, layer.N_y - 1))
        re = self.count_re(layer)
        for i in range(layer.N_y - 1):
            for j in range(layer.N_x - 1):
                if i == self.y_well and j == self.x_well:
                    wi[i][j] = -(layer.count_ro(pressure_matrix[i][j]) * 2 * math.pi * layer.k * layer.z / layer.mu) * (1 / math.log(re / self.r_w))
        return wi


class Solver:
    def __init__(self):
        self.delta_0 = 1000.0  # nachalnoe priblizhenie
        self.delta_max = 10 ** (-3)  # zadannaya tochnost
        self.delta_gr = 1000.0
        self.ksi_max = 10 ** (-3)
        self.tau = 86400.0  # s- time step
        # time_max = 365.25 * 4 * tau
        self.time_max = 365.25 * 1 * self.tau
        self.t_out_fields = 10 * self.tau  # frecuency of writing

    @staticmethod
    def count_tx(pressure_matrix, layer):
        tx_plus = np.zeros((layer.N_x - 1, layer.N_y - 1))
        tx_minus = np.zeros((layer.N_x - 1, layer.N_y - 1))
        for i in range(layer.N_y - 1):
            for j in range(layer.N_x - 1):
                if j == layer.N_x - 2:
                    tx_plus[i][j] = 0.0
                else:
                    tx_plus[i][j] = (layer.k_x * (layer.h * layer.z) / (layer.mu * layer.h)) * layer.count_ro(max(pressure_matrix[i][j], pressure_matrix[i][j + 1]))

            for j in range(layer.N_x - 1):
                if j == 0:
                    tx_minus[i][j] = 0.0
                else:
                    tx_minus[i][j] = (layer.k_x * (layer.h * layer.z) / (layer.mu * layer.h)) * layer.count_ro(max(pressure_matrix[i][j], pressure_matrix[i][j - 1]))
        tx_list = [tx_minus, tx_plus]
        return tx_list

    @staticmethod
    def count_ty(pressure_matrix, layer):
        ty_plus = np.zeros((layer.N_x - 1, layer.N_y - 1))
        ty_minus = np.zeros((layer.N_x - 1, layer.N_y - 1))
        for j in range(layer.N_x - 1):
            for i in range(layer.N_y - 1):
                if i == layer.N_y - 2:
                    ty_plus[i][j] = 0.0
                else:
                    ty_plus[i][j] = (layer.k_y * (layer.h * layer.z) / (layer.mu * layer.h)) * layer.count_ro(max(pressure_matrix[i][j], pressure_matrix[i + 1][j]))
            for i in range(layer.N_y - 1):
                if i == 0:
                    ty_minus[i][j] = 0.0
                else:
                    ty_minus[i][j] = (layer.k_y * (layer.h * layer.z) / (layer.mu * layer.h)) * layer.count_ro(max(pressure_matrix[i][j], pressure_matrix[i - 1][j]))
        ty_list = [ty_minus, ty_plus]
        return ty_list

    @staticmethod
    def check_diagonal(layer, a, b, c, g, f):
        for i in range(layer.N_y - 1):
            for j in range(layer.N_x - 1):
                if math.fabs(a[i][j]) < math.fabs(b[i][j]) + math.fabs(c[i][j]) + math.fabs(g[i][j]) + math.fabs(
                        f[i][j]):
                    print("Не сходится в ячейке  " + str(i) + " " + str(j))

    @staticmethod
    def count_r(layer, pressure_fictious, c, g, a, f, b, d):
        r = np.zeros((layer.N_x - 1, layer.N_y - 1))
        for i in range(layer.N_y - 1):
            r[i] = c[i] * pressure_fictious[i + 1][:-2] + g[i] * pressure_fictious[i][1:-1] + a[i] * pressure_fictious[i + 1][1:-1] + f[i] * pressure_fictious[i + 2][1:-1] + b[i] * pressure_fictious[i + 1][2:] - d[i]  # Check!!!!!?????
        return r

    @staticmethod
    def count_tolerance(matrix):
        matrix_abs = np.absolute(matrix)
        return np.amax(matrix_abs)

    @staticmethod
    #GZ method
    def count_ksi(layer, ksi, r, r_diff, r_x_diff_minus, r_x_diff_plus, r_y_diff_minus, r_y_diff_plus):
        for i in range(layer.N_y - 1):
            ksi[i + 1][1:-1] = (1 / r_diff[i]) * (-r[i] - r_x_diff_plus[i] * ksi[i + 1][2:] - r_y_diff_plus[i] * ksi[i + 2][1:-1] - r_x_diff_minus[i] * ksi[i + 1][:-2] - r_y_diff_minus[i] * ksi[i][1:-1])
        return ksi

    @staticmethod
    #Nonlinear GZ method
    #Going up and down
    def thomas_method_up_down(layer, ksi, r, r_diff, r_x_diff_minus, r_x_diff_plus, r_y_diff_minus, r_y_diff_plus, index):
        d_new = -r[index] - r_y_diff_minus[index] * ksi[index][1:-1] - r_y_diff_plus[index] * ksi[index + 2][1:-1]
        w_1 = r_diff[index][0]
        q_list = []
        w_list = [w_1]
        for j in range(1, len(r_diff[index])):
            q = r_x_diff_plus[index][j - 1] / w_list[j - 1]
            q_list.append(q)
            w = r_diff[index][j] - r_x_diff_minus[index][j] * q_list[j - 1]
            w_list.append(w)

        g_list = [d_new[0] / w_list[0]]

        for j in range(1, len(d_new)):
            g = (d_new[j] - r_x_diff_minus[index][j] * g_list[j - 1]) / w_list[j]
            g_list.append(g)

        list_ksi = [g_list[-1]]
        for j in range(len(g_list) - 2, -1, -1):
            ksi_j = g_list[j] - q_list[j] * list_ksi[-1]
            list_ksi.append(ksi_j)

        list_ksi.reverse()
        for j in range(layer.N_x - 1):
            ksi[index + 1][j + 1] = list_ksi[j]
        return ksi

    @staticmethod
    def thomas_method_left_right(layer, ksi, r, r_diff, r_x_diff_minus, r_x_diff_plus, r_y_diff_minus, r_y_diff_plus, index):
        d_new = -r[:, index] - r_x_diff_minus[:, index] * ksi[:, index][1:-1] - r_x_diff_plus[:, index] * ksi[:, index + 2][1:-1]
        w_1 = r_diff[:, index][0]
        q_list = []
        w_list = [w_1]
        for j in range(1, len(r_diff[:][index])):#!!!!!!
            q = r_y_diff_plus[j - 1][index] / w_list[j - 1]
            q_list.append(q)
            w = r_diff[j][index] - r_y_diff_minus[j][index] * q_list[j - 1]
            w_list.append(w)

        g_list = [d_new[0] / w_list[0]]

        for j in range(1, len(d_new)):
            g = (d_new[j] - r_y_diff_minus[j][index] * g_list[j - 1]) / w_list[j]
            g_list.append(g)

        list_ksi = [g_list[-1]]
        for j in range(len(g_list) - 2, -1, -1):
            ksi_j = g_list[j] - q_list[j] * list_ksi[-1]
            list_ksi.append(ksi_j)

        list_ksi.reverse()
        for j in range(layer.N_y - 1):
            ksi[j + 1][index + 1] = list_ksi[j]
        return ksi


if __name__ == "__main__":
    layer = Layer()
    solver = Solver()
    well = Well()
    N_x = layer.N_x
    N_y = layer.N_y
    p = np.ones((N_x - 1, N_y - 1)) * layer.P_init
    p_gr = np.ones((N_y + 1, N_x + 1)) * layer.P_init
    delta = np.ones((N_y + 1, N_x + 1)) * solver.delta_0
    delta_k = np.zeros((N_y + 1, N_x + 1))
    ksi = delta.copy()
    ksi_k = np.ones((N_y + 1, N_x + 1))
    q_list = []
    time_list = []
    time = solver.tau
    counter = 1

    while time < solver.time_max:
        delta = np.ones((N_y + 1, N_x + 1)) * solver.delta_0
        delta_k = np.zeros((N_y + 1, N_x + 1))
        p_n = p.copy()
        while solver.count_tolerance(delta_k - delta) > solver.delta_max:
            ksi = np.ones((N_y + 1, N_x + 1)) * solver.delta_0
            beta = layer.count_beta(p)
            wi = well.count_wi(p, layer)
            c, b = solver.count_tx(p, layer)
            g, f = solver.count_ty(p, layer)
            a = -(b + c + g + f) + wi - layer.V_ij * beta / solver.tau
            d = wi * well.pressure_w - p_n * (layer.V_ij * beta / solver.tau)

            solver.check_diagonal(layer, a, b, c, g, f)
            r_x_diff_minus = c.copy()
            r_x_diff_plus = b.copy()
            r_y_diff_minus = g.copy()
            r_y_diff_plus = f.copy()
            r_diff = a.copy()
            r = solver.count_r(layer, p_gr, c, g, a, f, b, d)
            delta_k = delta.copy()

            while solver.count_tolerance(ksi - ksi_k) > solver.delta_max:
                ksi_k = ksi.copy()
                for i in range(layer.N_y - 1):
                    ksi = solver.thomas_method_up_down(layer, ksi, r, r_diff, r_x_diff_minus, r_x_diff_plus, r_y_diff_minus, r_y_diff_plus, i)
                for i in range(layer.N_y - 2, -1, -1):
                    ksi = solver.thomas_method_up_down(layer, ksi, r, r_diff, r_x_diff_minus, r_x_diff_plus, r_y_diff_minus,r_y_diff_plus, i)

                for i in range(layer.N_x - 1):
                    ksi = solver.thomas_method_left_right(layer, ksi, r, r_diff, r_x_diff_minus, r_x_diff_plus, r_y_diff_minus,
                                                       r_y_diff_plus, i)

                for i in range(layer.N_x - 2, -1, -1):
                    ksi = solver.thomas_method_left_right(layer, ksi, r, r_diff, r_x_diff_minus, r_x_diff_plus, r_y_diff_minus,
                                                       r_y_diff_plus, i)

                tolerance = solver.count_tolerance(ksi - ksi_k)
                #print(tolerance)

            delta = ksi.copy()
            tolerance = solver.count_tolerance(delta - delta_k)
            print(str(tolerance) + " meow")

            for i in range(N_y - 1):
                for j in range(N_x - 1):
                    p[i][j] += delta[i + 1][j + 1]

            p_gr += delta

        if counter % 1 == 0:
            name = 'graph_' + str(counter)
            fig = plt.figure()
            ax = fig.gca(projection='3d')

            # Make data.
            X = np.arange(layer.x_0 + layer.h / 2, layer.x_N, layer.h)
            Y = np.arange(layer.y_0 + layer.h / 2, layer.y_N, layer.h)
            X, Y = np.meshgrid(X, Y)

            # Plot the surface.
            surf = ax.plot_surface(X, Y, p, cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False)
            plt.savefig(name)
            plt.clf()

            x = np.arange(layer.x_0 + layer.h / 2, layer.x_N, layer.h)
            file = open('pressure_10-days.txt', 'w')
            for j in range(len(p[0])):
                file.write(str(x[j]) + ' ' + str(p[well.y_well][j]) + '\n')
            file.close()

        # q = wi[45][45] * (p[45][45] - well.pressure_w)
        # q_list.append(q)
        # time_list.append(time)
        #
        # print(time / 86400.0)
        time += solver.tau
        counter += 1

    # file = open('q_time.txt', 'w')
    # for i in range(len(q_list)):
    #     file.write(str(time_list[i]) + '  ' + str(q_list[i]) + '\n')
    # file.close()
    #
    # plt.plot(time_list, q_list)
    # plt.show()
    #
    # #
    #
    # stuff to run always here such as class/def

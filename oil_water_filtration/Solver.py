import numpy as np


class SolverSlau():
    def __init__(self, equation_count): #equation_count - количество уравнений равно количеству неизвестных
        self.e_count = equation_count
        self.coefficient_matrix = np.matrix((self.e_count, self.e_count))
        self.result_vector = np.matrix(self.e_count)
        self.nevyaz_vector = np.matrix(self.e_count)

    def add_coefficient(self, i, j, coefficient):
        self.coefficient_matrix[i, j] += coefficient

    def set_coefficient(self, i, j, value):
        self.coefficient_matrix[i, j] = value

    def set_zero(self):
        self.coefficient_matrix * 0.0

    def add_nevyaz(self, i, coefficient):
        self.nevyaz_vector[i] += coefficient


    def solve_thomas_method(self, list_d, P0, PN):
        w_1 = self.coefficient_matrix[0, 0]
        q_list = []
        w_list = [w_1]
        for j in range(1, self.e_count):
            q = self.coefficient_matrix[j - 1, j] / w_list[j - 1]
            q_list.append(q)
            w = self.coefficient_matrix[j, j] - self.coefficient_matrix[j, j - 1] * q_list[j - 1]
            w_list.append(w)

        g_list = [list_d[0] / w_list[0]]
        for j in range(1, self.e_count):
            g = (list_d[j] - self.coefficient_matrix[j, j - 1] * g_list[j - 1]) / w_list[j]
            g_list.append(g)

        #list_P = [g_list[-1]]
        self.result_vector[-1] = g_list[-1]
        for j in range(len(g_list) - 2, -1, -1):
            res = g_list[j] - q_list[j] * self.result_vector[-1]
            self.result_vector.append(res)


        # Add usloviya na granitsah
        #list_P.append(P0)
        #list_P.reverse()
        #list_P.append(PN)

    def get_result(self):
        return self.result_vector



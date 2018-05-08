import numpy as np


class SolverSlau():
    def __init__(self, equation_count): #equation_count - количество уравнений равно количеству неизвестных
        self.e_count = equation_count
        self.coefficient_matrix = np.zeros((self.e_count, self.e_count), dtype=object)
        #self.result_vector = np.array(self.e_count)
        self.result_vector = []
        self.nevyaz_vector = np.zeros(self.e_count)

    def add_coefficient(self, i, j, coefficient):
        if (i == 0):
            print(i, j, coefficient )
        self.coefficient_matrix[i, j] += coefficient

    def set_coefficient(self, i, j, value):
        self.coefficient_matrix[i, j] = value

    def set_zero(self):
        self.coefficient_matrix = np.zeros((self.e_count, self.e_count))
        self.nevyaz_vector = np.zeros(self.e_count)

    def add_nevyaz(self, i, coefficient):
        self.nevyaz_vector[i] += coefficient

    def set_matrix_coefficients(self, i, j, pressure_oil_i, pressure_oil_j, coefficient):
        if i < self.e_count and i >= 0:
            self.add_coefficient(i, i, -coefficient)
            # self.add_nevyaz(i, coefficient * pressure_oil[i + 1])
            # self.add_nevyaz(i, -coefficient * pressure_oil[j + 1])
            self.add_nevyaz(i, coefficient * pressure_oil_i)
            self.add_nevyaz(i, -coefficient * pressure_oil_j)

            if j < self.e_count and j >= 0:
                self.add_coefficient(i, j, coefficient)

    def set_matrix_coefficients_ss(self, i, j, coefficient):
        if i < self.e_count and i >= 0:
            self.add_coefficient(i, i, -coefficient)
            if j < self.e_count and j >= 0:
                self.add_coefficient(i, j, coefficient)



    #for a and d
    def set_matrix_coeffitients_a_d(self, i, j, pressure_oil, pressure_oil_n, coefficient):
        #A_c1_p
        if i < self.e_count and i >= 0:
            self.add_coefficient(i, i, -coefficient)
            self.add_nevyaz(i, coefficient * pressure_oil)
            self.add_nevyaz(i, -coefficient * pressure_oil_n)

    def solve_thomas_method(self, left_gr, right_gr):
        w_1 = self.coefficient_matrix[0, 0]
        q_list = []
        w_list = [w_1]
        for j in range(1, self.e_count):
            q = self.coefficient_matrix[j - 1, j] / w_list[j - 1]
            q_list.append(q)
            w = self.coefficient_matrix[j, j] - self.coefficient_matrix[j, j - 1] * q_list[j - 1]
            w_list.append(w)

        g_list = [self.nevyaz_vector[0] / w_list[0]]
        for j in range(1, self.e_count):
            g = (self.nevyaz_vector[j] - self.coefficient_matrix[j, j - 1] * g_list[j - 1]) / w_list[j]
            g_list.append(g)

        self.result_vector.append(g_list[-1])
        for j in range(len(g_list) - 2, -1, -1):
            res = g_list[j] - q_list[j] * self.result_vector[-1]
            self.result_vector.append(res)


        # Add usloviya na granitsah
        self.result_vector.append(left_gr)
        self.result_vector.reverse()
        self.result_vector.append(right_gr)

    def solve_thomas_matrix_method(self, left_gr_water, right_gr_water, left_gr_oil, right_gr_oil):
        w_1 = self.coefficient_matrix[0, 0]
        #q_1 = np.linalg.inv(w_1) * self.coefficient_matrix[0, 1]
        w_list = [w_1]
        q_list = []

        for i in range(1, self.e_count):
            q_i = np.linalg.inv(w_list[-1]) * self.coefficient_matrix[i - 1, i]
            q_list.append(q_i)
            w_i = self.coefficient_matrix[i, i] - self.coefficient_matrix[i, i - 1] * q_list[i - 1]
            w_list.append(w_i)


        g_1 = np.matmul(np.linalg.inv(w_list[0]),self.nevyaz_vector[0])
        g_list = [g_1]
        for i in range(1, self.e_count):
            b = self.nevyaz_vector[i] - np.matmul(self.coefficient_matrix[i, i - 1] , g_list[i - 1])
            g_i = np.matmul(np.linalg.inv(w_list[i]), b)
            g_list.append(g_i)

        self.result_vector.append(g_list[-1])
        for i in range(len(g_list) - 2, -1, -1):
            u_i = g_list[i] - np.matmul(q_list[i], self.result_vector[-1])
            self.result_vector.append(u_i)

        #self.result_vector.append([left_gr_water, left_gr_oil])
        self.result_vector.reverse()
        #self.result_vector.append([right_gr_water, right_gr_oil])


    def get_result(self):
        return self.result_vector

    def clear_result(self):
        self.result_vector = []

#
# a_list = [-3.8300546467228913e-08, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09, -6.963735738136903e-09]
#
# b_list = [3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09]
#
# c_list = [3.4818678587880006e-08, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09, 3.4818678587880006e-09]
#
# d_list = [-0.17640013039584707, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20, -8.772599362945779e-20]
#
# delta_k = [0.0, 5061085.527884434, 5009440.836615181, 4957796.174927355, 4906151.5425159875, 4854506.939076111, 4802862.364302755, 4751217.817890952, 4699573.2995357355, 4647928.808932137, 4596284.345775187, 4544639.90975992, 4492995.500581369, 4441351.117934564, 4389706.76151454, 4338062.431016331, 4286418.126134968, 4234773.846565486, 4183129.592002917, 4131485.362142296, 4079841.1566786566, 4028196.975307032, 3976552.8177224565, 3924908.6836199644, 3873264.572694591, 3821620.48464137, 3769976.4191553355, 3718332.375931523, 3666688.3546649665, 3615044.355050702, 3563400.376783765, 3511756.41955919, 3460112.4830720127, 3408468.5670172684, 3356824.6710899933, 3305180.794985223, 3253536.9383979933, 3201893.10102334, 3150249.2825563, 3098605.482691909, 3046961.7011252036, 2995317.93755122, 2943674.1916649947, 2892030.463161564, 2840386.7517359657, 2788743.057083236, 2737099.3788984125, 2685455.7168765315, 2633812.0707126306, 2582168.440101747, 2530524.8247389183, 2478881.224319182, 2427237.638537575, 2375594.0670891353, 2323950.5096689006, 2272306.9659719085, 2220663.4356931974, 2169019.918527805, 2117376.414170769, 2065732.9223171277, 2014089.4426619194, 1962445.9749001823, 1910802.518726955, 1859159.0738372754, 1807515.6399261819, 1755872.2166887133, 1704228.803819908, 1652585.4010148048, 1600942.0079684425, 1549298.62437586, 1497655.2499320959, 1446011.8843321889, 1394368.527271178, 1342725.1784441026, 1291081.8375460012, 1239438.5042719129, 1187795.1783168768, 1136151.8593759325, 1084508.547144119, 1032865.2413164752, 981221.9415880408, 929578.6476538549, 877935.359208957, 826292.0759483862, 774648.7975671823, 723005.5237603845, 671362.2542230324, 619718.9886501654, 568075.7267368229, 516432.46817804454, 464789.2126688699, 413145.9599043385, 361502.70957949007, 309859.46138936415, 258216.2150290003, 206572.9701934382, 154929.72657771746, 103286.48387687776, 51643.24178595872, 0.0]
#
# solverSlau = SolverSlau(len(a_list))
# #a = np.matrix((98, 98))
# for i in range(solverSlau.e_count):
#     solverSlau.coefficient_matrix[i, i] = a_list[i]
#     solverSlau.nevyaz_vector[i] = d_list[i]
#
# for i in range(solverSlau.e_count - 1):
#     solverSlau.coefficient_matrix[i, i + 1] = b_list[i]
#
# for i in range(1, solverSlau.e_count):
#     solverSlau.coefficient_matrix[i, i - 1] = c_list[i]
#
# solverSlau.solve_thomas_method(0.0, 0.0)
# print("meow")
#

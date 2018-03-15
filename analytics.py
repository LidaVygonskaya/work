from matplotlib import pyplot as plt
import math
import numpy as np
# coding: utf8


class Integrals:
    def __init__(self):
        self.atm = 101325.0
        self.k = (9.868233 * (10 ** (-13))) * 10 ** (-3)
        self.fi = 0.2
        self.mu = 10 ** (-3)  # Pa * s
        self.c_f = (10 ** (-4)) / self.atm
        self.z = 5.05
        self.hi = self.k / (self.fi * self.mu * self.c_f)
        self.r_w = 0.108  # meters
        self.r_w_max = 270.0
        self.p_w = 30 * self.atm
        self.p_0 = 80 * self.atm
        self.N = 1000
        self.h = (self.r_w_max - self.r_w) / self.N

    #r - радиус скважины
    def f(self, r, t):
        return (1 /(r)) * math.exp(-r ** 2 / (4.0 * self.hi * t))

    def count_f_integral(self, r, t):
        step = (self.r_w_max - r) / self.N
        x_setka = np.arange(r, self.r_w_max, step)
        x_setka = np.append(x_setka, [self.r_w_max])

        summ = 0.0
        integral = 0.0
        for i in range(self.N):
            summ += 0.5 * (self.f(x_setka[i], t) + self.f(x_setka[i + 1], t)) * (x_setka[i + 1] - x_setka[i])
        integral += summ
        return integral

    def count_q(self, t):
        return -(2.0 * math.pi * self.k * self.z * 1000.0 / self.mu) * ((self.p_w - self.p_0) / self.count_f_integral(self.r_w, t)) * math.exp(-(self.r_w ** 2.0 / (4.0 * self.hi * t)))

    def count_pressure(self, r, t):
        return self.p_0 + ((self.p_w - self.p_0) / self.count_f_integral(self.r_w, t)) * self.count_f_integral(r, t)


#Аналитика для q
#Численное решение
q_list = []
time_list = []
file = open("q_time.txt", 'r')
for line in file.readlines():
    line = line.rstrip()
    s = line.split('  ')
    time_list.append(float(s[0]) / 86400.0)
    q_list.append(float(s[1]) * (-1))

#Аналитическое решение
integ = Integrals()
q_list_anal = []
for time in time_list:
    q = integ.count_q(time)
    q_list_anal.append(q)

plt.plot(time_list, q_list)
plt.plot(time_list, q_list_anal)
plt.xlabel('time, days')
plt.ylabel('Q')
plt.show()

file = open('q_time_anal.txt', 'w')
for i in range(len(q_list_anal)):
    file.write(str(time_list[i]) + '  ' + str(q_list_anal[i]) + '\n')
file.close()

#Аналитика для P
#Численное решение
p_list = []
x_list = []
file = open("pressure_10-days.txt", 'r')
for line in file.readlines():
    line = line.rstrip()
    s = line.split(' ') #разделитель пробел

    x_list.append(float(s[0]))
    p_list.append(float(s[1]))

#Аналитическое решение
integ = Integrals()

x = np.arange(0.108, 250, (250 - 0.108) / 1000)
x = np.append(x, [250])
p_list_anal = []

for elem in x[:-1]:
    p = integ.count_pressure(elem, time_list[9] * 86400.0)
    p_list_anal.append(p)


plt.plot(x_list, p_list)
plt.plot(x[:-1] + 230, p_list_anal)
plt.xlabel('x, meters')
plt.ylabel('P, Pa')
plt.show()

file = open('pressure-10-days_anal.txt', 'w')
for i in range(len(q_list_anal)):
    file.write(str(x[:-1][i]) + '  ' + str(p_list_anal[i]) + '\n')
file.close()


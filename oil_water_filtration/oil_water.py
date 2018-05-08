import math
import matplotlib.pyplot as plt
import numpy as np

from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
from oil_water_filtration.OilWaterIMPES import OilWaterIMPES
from oil_water_filtration.Solver import SolverSlau


layer = Layer()
oil_water_impes = OilWaterIMPES()
solver_slau = SolverSlau(layer.N - 2)

cell_container = CellContainer(layer.N, layer)
cell_container.initialize_cells()

flow_array = []
for i in range(cell_container.get_len() - 1):
    cells = cell_container.get_cells()
    flow_array.append(Flow(cells[i], cells[i + 1]))

time = oil_water_impes.tau
counter = 1

while time < oil_water_impes.time_max:

    delta_list_k = [0.0] + [oil_water_impes.delta_0 for i in range(layer.N - 2)] + [0.0]

    for cell in cell_container.get_cells():
        cell.get_cell_state_n_plus().set_pressure_cap(0.0)
        cell.get_cell_state_n().set_equals_to(cell.get_cell_state_n_plus())

    while oil_water_impes.count_norm(delta_list_k) > oil_water_impes.delta_max:
        solver_slau.set_zero()
        oil_water_impes.recount_properties(cell_container)
        oil_water_impes.count_flows(flow_array)
        oil_water_impes.generate_matrix(flow_array, cell_container, solver_slau, oil_water_impes)

        solver_slau.solve_thomas_method(0.0, 0.0)
        delta_list_k = solver_slau.get_result()
        solver_slau.clear_result()
        oil_water_impes.update_pressure(cell_container, delta_list_k)

    oil_water_impes.recount_properties(cell_container)
    oil_water_impes.count_flows(flow_array)
    oil_water_impes.update_saturation(cell_container, flow_array, oil_water_impes)
    #cells = cell_container.get_cells()
    #cells[-2].get_state_n = cells[-1].

    if oil_water_impes.count_norm([cell.get_cell_state_n_plus().get_s_water() - cell.get_cell_state_n().get_s_water() for cell in cell_container.get_cells()]) > 0.01:
        print("S_n+1 - S_n > 0.1")
        oil_water_impes.tau = oil_water_impes.tau / 2.0
        #the same
        for cell in cell_container.get_cells():
            cell.get_cell_state_n_plus().set_equals_to(cell.get_cell_state_n())

    else:
        print(time)
        time += oil_water_impes.tau
        oil_water_impes.tau = max(oil_water_impes.tau * 2.0, oil_water_impes.tau_default)
        if 125 * 86400.0 - time - oil_water_impes.tau <= 0:
            oil_water_impes.tau = 125 * 86400.0 - time #- solver.tau

        if abs(time - 86400.0 * 125) < 1.e-16:
            x = np.arange(layer.x_0, layer.x_N, layer.h)
            x = np.append(x, layer.x_N)
            #
            #plt.plot(x, P_water, 'bo', x, P_water, 'k')
            #plt.title(str(time / 86400.0) + 'days')
            #plt.xlabel('x')
            #plt.ylabel('P_water')
            #plt.show()
            # plt.plot(x, s_water, 'bo', x, s_water, 'k')
            # plt.title(str(time / 86400.0) + 'days')
            # plt.xlabel('x')
            # plt.ylabel('S_water')
            # plt.show()
            file = open('s(x)_count_125_days.txt', 'w')
            for i in range(cell_container.get_len()):
                state = cell_container.get_cells()[i].get_cell_state_n_plus()
                file.write(str(x[i]) + ' ' + str(state.get_s_water()) + '\n')

            file.close()
        counter += 1


























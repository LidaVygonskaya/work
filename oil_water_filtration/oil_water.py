from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
from oil_water_filtration.OilWaterIMPES import OilWaterIMPES
from oil_water_filtration.Solver import SolverSlau
import numpy as np
from matplotlib import pyplot as plt

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

    delta_list_k = oil_water_impes.generate_delta_k(layer)

    for cell in cell_container.get_cells():
        cell.get_cell_state_n().set_equals_to(cell.get_cell_state_n_plus())

    while oil_water_impes.count_norm(delta_list_k) > oil_water_impes.delta_max:
        solver_slau.set_zero()
        oil_water_impes.recount_properties(cell_container)
        oil_water_impes.count_flows(flow_array)
        oil_water_impes.generate_matrix(flow_array, cell_container, solver_slau)

        oil_water_impes.solve_slau(solver_slau)
        delta_list_k = solver_slau.get_result()
        solver_slau.clear_result()
        #oil_water_impes.update_pressure(cell_container, delta_list_k)

        if oil_water_impes.check_pressure_convergence(cell_container, delta_list_k):
            oil_water_impes.tau = oil_water_impes.tau / 2.0
            for cell in cell_container.get_cells():
                cell.get_cell_state_n_plus().set_equals_to(cell.get_cell_state_n())
        else:
            oil_water_impes.update_pressure(cell_container, delta_list_k)

    oil_water_impes.recount_properties(cell_container)
    oil_water_impes.count_flows(flow_array)
    oil_water_impes.update_saturation(cell_container, flow_array)
    oil_water_impes.update_pressure_cap(cell_container)

    # x = np.arange(layer.x_0, layer.x_N, layer.h)
    # x = np.append(x, layer.x_N)
    # s_wat = []
    # s_oil = []
    #
    # for cell in cell_container.get_cells():
    #     s_wat.append(cell.get_cell_state_n_plus().get_s_water())
    #
    # plt.plot(x, s_wat, 'r')
    # plt.show()

    print("CURRENT TIME ", time / 86400)
    if time / 86400 > 1.0302:
        print("meme")

    if oil_water_impes.check_saturation_convergence(cell_container):
        oil_water_impes.tau = oil_water_impes.tau / 2.0
        for cell in cell_container.get_cells():
            cell.get_cell_state_n_plus().set_equals_to(cell.get_cell_state_n())
        if time / 86400 > 16.15:
            x = np.arange(layer.x_0, layer.x_N, layer.h)
            x = np.append(x, layer.x_N)
            s_wat = []
            s_oil = []
            for cell in cell_container.get_cells():
                s_wat.append(cell.get_cell_state_n_plus().get_s_water())
                s_oil.append(cell.get_cell_state_n_plus().get_s_oil())

            #plt.plot(x, s_wat, 'r')
            #plt.plot(x, s_oil)
            #plt.plot(x[:-1], s_oil_check)
            #plt.show()

    else:
        # t_days = 50
        # print(time / 86400.0, '        ', oil_water_impes.tau)
        # oil_water_impes.tau = min(oil_water_impes.tau * 2.0, oil_water_impes.tau_default)
        # if (time / 86400 == t_days) and (time % 86400 == 0):
        #     oil_water_impes.show_results(time, layer, cell_container)
        # time += oil_water_impes.tau
        # counter += 1


        t_days = 50
        print(time / 86400.0, '        ', oil_water_impes.tau)
        oil_water_impes.tau = min(oil_water_impes.tau * 2.0, oil_water_impes.tau_default)
        if t_days * 86400.0 - time - oil_water_impes.tau <= 0:
            oil_water_impes.tau = t_days * 86400.0 - time

        if abs(time - 86400.0 * t_days) < 1.e-16:
            oil_water_impes.show_results(time, layer, cell_container)
        time += oil_water_impes.tau
        counter += 1
        #
























from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
from oil_water_filtration.OilWaterSS import OilWaterSS
from oil_water_filtration.Solver import SolverSlau
from matplotlib import pyplot as plt
import numpy as np

layer = Layer()
oil_water_ss = OilWaterSS()
solver_slau = SolverSlau(layer.N - 2)

cell_container = CellContainer(layer.N, layer)
cell_container.initialize_cells()

flow_array = []
for i in range(cell_container.get_len() - 1):
    cells = cell_container.get_cells()
    flow_array.append(Flow(cells[i], cells[i + 1]))

time = oil_water_ss.tau
counter = 1

while time < oil_water_ss.time_max:
    if oil_water_ss.tau < 86400.0:
        oil_water_ss.tau = oil_water_ss.tau * 2.0

    delta_list_k = oil_water_ss.generate_delta_k(layer)
    for cell in cell_container.get_cells():
        cell.get_cell_state_n().set_equals_to(cell.get_cell_state_n_plus())

    for cell in cell_container.get_cells()[1:-1]:
        cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n().get_pressure_water() + 1000.0)
        cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n().get_pressure_oil() + 1000.0)

    while oil_water_ss.count_norm(delta_list_k) > oil_water_ss.delta_max:
        print(oil_water_ss.count_norm(delta_list_k))
        solver_slau.set_zero()
        oil_water_ss.reshape_matrices(solver_slau)
        oil_water_ss.recount_properties(cell_container)
        oil_water_ss.count_flows(flow_array)
        oil_water_ss.generate_matrix(flow_array, cell_container, solver_slau)
        oil_water_ss.solve_slau(solver_slau)
        delta_list_k = solver_slau.get_result()
        solver_slau.clear_result()

        #Проверка на сходимость
        #if oil_water_ss.check_pressure_convergence(cell_container, delta_list_k):
        if oil_water_ss.check_pressure_convergence(cell_container, delta_list_k):
            oil_water_ss.tau = oil_water_ss.tau / 2.0
            for cell in cell_container.get_cells():
                cell.get_cell_state_n_plus().set_equals_to(cell.get_cell_state_n())
        else:
            oil_water_ss.update_pressure(cell_container, delta_list_k)


        #Improvements
        #or cell in cell_container.get_cells()[1:-1]:
        #   state_n_plus = cell.get_cell_state_n_plus()
        #   state_n = cell.get_cell_state_n()
        #   cell_layer = cell.layer

        #   # Recount pressure_cap
        #   p_oil = state_n_plus.get_pressure_oil()
        #   p_water = state_n_plus.get_pressure_water()
        #   state_n_plus.set_pressure_cap(p_oil - p_water)

        #   #Recount s_water
        #   s_water_graph = layer.count_s_water_graph()
        #   s_der = cell_layer.count_s_water_graph_der_from_p(s_water_graph, state_n.get_pressure_cap())
        #   s_wat_new = state_n.get_s_water() + s_der * (state_n_plus.get_pressure_cap() - state_n.get_pressure_cap())
        #   state_n_plus.set_s_water(s_wat_new)

        #   #Recount p_cap again
        #   state_n_plus.set_pressure_cap(cell_layer.count_pressure_cap(s_water_graph, state_n_plus.get_s_water()))

        #   #Update Pressure_water
        #   state_n_plus.set_pressure_water(state_n_plus.get_pressure_oil() - state_n_plus.get_pressure_cap())

    #oil_water_ss.generate_matrix(flow_array, cell_container, solver_slau)

    oil_water_ss.recount_properties(cell_container)
    oil_water_ss.update_saturation(cell_container, flow_array)

    #RES
    # x = np.arange(layer.x_0, layer.x_N, layer.h)
    # x = np.append(x, layer.x_N)
    # s_wat = []
    # for cell in cell_container.get_cells():
    #     s_wat.append(cell.get_cell_state_n_plus().get_s_water())
    #
    # plt.plot(x, s_wat, 'r')
    # plt.show()

    print("COUNTER " + str(counter))

    if counter == 125 or counter == 364 or counter == 50 or counter == 1 or counter == 40:
        oil_water_ss.show_results(counter, layer, cell_container)

    counter += 1
    print('DAYS: ', time / 86400, '        ','STEP: ',  oil_water_ss.tau)
    time += oil_water_ss.tau


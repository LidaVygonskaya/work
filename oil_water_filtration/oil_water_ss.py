from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
from oil_water_filtration.OilWaterSS import OilWaterSS
from oil_water_filtration.Solver import SolverSlau
import numpy as np
import scipy.sparse.linalg as sc
from oil_water_filtration.Write import write_nevyaz,write
from matplotlib import pyplot as plt

layer = Layer()
oil_water_ss = OilWaterSS()
solver_slau = SolverSlau(layer.N - 2)
solver_slau.coefficient_matrix = np.zeros((98, 98, 2, 2), float)
solver_slau.nevyaz_vector = np.zeros((98, 2, 1))

cell_container = CellContainer(layer.N, layer)
cell_container.initialize_cells()


for cell in cell_container.get_cells()[1:]:
    cell.get_cell_state_n().set_s_water(0.1)
    #cell.get_cell_state_n_plus().set_s_water(10 ** (-4))

flow_array = []
for i in range(cell_container.get_len() - 1):
    cells = cell_container.get_cells()
    flow_array.append(Flow(cells[i], cells[i + 1]))


time = oil_water_ss.tau
counter = 1

while time < oil_water_ss.time_max:
    delta_list_k = [[0.0, 0.0]] + [[oil_water_ss.delta_0, oil_water_ss.delta_0] for i in range(layer.N - 2)] + [[0.0, 0.0]]
    for cell in cell_container.get_cells():
        cell.get_cell_state_n().set_equals_to(cell.get_cell_state_n_plus())

    # for cell in cell_container.get_cells()[1:-1]:
    #     cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n_plus().get_pressure_water() + 1000.0)
    #     cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n_plus().get_pressure_oil() + 1000.0)

    if counter == 1:
        for cell in cell_container.get_cells()[1:]:
            cell.get_cell_state_n().set_s_water(0.1)
            cell.get_cell_state_n().set_s_oil(1.0 - 0.1)
            cell.get_cell_state_n().set_s_oil(0.9 ** 2)
            cell.get_cell_state_n().set_s_water(0.1 ** 2)
            cell.get_cell_state_n_plus().set_s_water(0.1)
            cell.get_cell_state_n_plus().set_s_oil(0.9)




#delta_list_k_bicg = [0.0] + [oil_water_ss.delta_0 for i in range(layer.N - 2)] + [0.0]
    while oil_water_ss.count_matrix_norm(delta_list_k) > oil_water_ss.delta_max:
        print(oil_water_ss.count_matrix_norm(delta_list_k))
        solver_slau.set_zero()
        solver_slau.coefficient_matrix = np.zeros((98, 98, 2, 2), float)
        solver_slau.nevyaz_vector = np.zeros((98, 2, 1))
        oil_water_ss.recount_properties(cell_container)
        oil_water_ss.count_flows(flow_array)

        #generate matrix
        i = 0
        for flow in flow_array:
            oil_water_ss.count_b(solver_slau, flow)
            oil_water_ss.count_c(solver_slau, flow)
            i+=1

        coeff_mat = write(solver_slau.coefficient_matrix)
        nevyaz_mat = write_nevyaz(solver_slau.nevyaz_vector)

        for cell in cell_container.get_cells():
            oil_water_ss.count_a_ss(solver_slau, cell)

        coeff_mat = write(solver_slau.coefficient_matrix)
        nevyaz_mat = write_nevyaz(solver_slau.nevyaz_vector)
        #generate_matrix

        #Улучшение сходимости
        # for i in range(coeff_mat.shape[0]):
        #     a = coeff_mat[i, i]
        #     nevyaz_mat[i] = nevyaz_mat[i] / a
        #     for j in range(coeff_mat.shape[1]):
        #         coeff_mat[i, j] = coeff_mat[i, j] / a

        #delta_list_k_bicg = sc.spsolve(coeff_mat, nevyaz_mat)
        #delta_list_k_bicg = delta_list_k_bicg[0]
        solver_slau.solve_thomas_matrix_method(0.0, 0.0, 0.0, 0.0)
        delta_list_k = solver_slau.get_result()

        solver_slau.clear_result()

        #update pressure
        for i in range(1, cell_container.get_len() - 1):
            cell = cell_container.get_cells()[i]
            cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n_plus().get_pressure_water() + delta_list_k[i - 1][0][0])
            cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n_plus().get_pressure_oil() + delta_list_k[i - 1][1][0])

        p_wat = []
        p_oil = []
        for cell in cell_container.get_cells():
            p_wat.append(cell.get_cell_state_n_plus().get_pressure_water())
            p_oil.append(cell.get_cell_state_n_plus().get_pressure_oil())

        x = np.arange(layer.x_0, layer.x_N, layer.h)
        x = np.append(x, layer.x_N)
        # plt.plot(x ,p_wat)
        # plt.plot(x, p_oil)
        # plt.show()


    #update everything
    # m = 1
    # for i in range(0, len(delta_list_k_bicg), 2):
    #     cell = cell_container.get_cells()[m]
    #     cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n_plus().get_pressure_water() + delta_list_k_bicg[i])
    #     m += 1
    #
    # m = 1
    # for i in range(1, len(delta_list_k_bicg), 2):
    #     cell = cell_container.get_cells()[m]
    #     cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n_plus().get_pressure_oil() + delta_list_k_bicg[i])
    #     m += 1

    p_wat = []
    p_oil = []
    x = np.arange(layer.x_0, layer.x_N, layer.h)
    x = np.append(x, layer.x_N)

    #update cells pressure_cap and saturation
    for cell in cell_container.get_cells()[1:-1]:
        p_wat.append(cell.get_cell_state_n_plus().get_pressure_water())
        p_oil.append(cell.get_cell_state_n_plus().get_pressure_oil())
        state_n_plus = cell.get_cell_state_n_plus()
        state_n_plus.set_pressure_cap(state_n_plus.get_pressure_oil() - state_n_plus.get_pressure_water())
        s_water_graph = cell.layer.count_s_water_graph()
        print(state_n_plus.get_pressure_cap())
        s_w = cell.layer.count_s_water(s_water_graph, state_n_plus.get_pressure_cap())
        cell.get_cell_state_n_plus().set_s_water(s_w)
        cell.get_cell_state_n_plus().set_s_oil(1.0 - s_w)

    print("COUNTER " + str(counter))
    if counter == 364:
        file = open('s(x)_count_364_days_ss_method.txt', 'w')
        for i in range(cell_container.get_len()):
            state = cell_container.get_cells()[i].get_cell_state_n_plus()
            file.write(str(x[i]) + ' ' + str(state.get_s_water()) + '\n')
        file.close()
    #plt.plot(x, p_wat)
    #plt.plot(x, p_oil)
    #plt.show()
    counter += 1
    time += oil_water_ss.tau


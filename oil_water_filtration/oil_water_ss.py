from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
from oil_water_filtration.OilWaterSS import OilWaterSS
from oil_water_filtration.Solver import SolverSlau
import numpy as np
import scipy.sparse.linalg as sc
from oil_water_filtration.Write import write_nevyaz,write

layer = Layer()
oil_water_ss = OilWaterSS()
solver_slau = SolverSlau(layer.N - 2)
solver_slau.coefficient_matrix = np.zeros((98, 98, 2, 2), float)
solver_slau.nevyaz_vector = np.zeros((98, 2, 1))

cell_container = CellContainer(layer.N, layer)
cell_container.initialize_cells()
#
# for cell in cell_container.get_cells():
#     cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n_plus().get_pressure_oil() + 1000.0)
#     cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n_plus().get_pressure_water() + 1000.0)
#
# for cell in cell_container.get_cells()[1:-1]:
#     cell.get_cell_state_n_plus().set_s_water(0.1)
#     cell.get_cell_state_n_plus().set_s_oil(0.9)


flow_array = []
for i in range(cell_container.get_len() - 1):
    cells = cell_container.get_cells()
    flow_array.append(Flow(cells[i], cells[i + 1]))

delta_list_k = [0.0] + [oil_water_ss.delta_0 for i in range(layer.N - 2)] + [0.0]
while oil_water_ss.count_norm(delta_list_k) > oil_water_ss.delta_max:
    print(oil_water_ss.count_norm(delta_list_k))
    solver_slau.set_zero()
    solver_slau.coefficient_matrix = np.zeros((98, 98, 2, 2), float)
    solver_slau.nevyaz_vector = np.zeros((98, 2, 1))
    oil_water_ss.recount_properties(cell_container)
    oil_water_ss.count_flows(flow_array)
    #generate matrix
    for flow in flow_array:
        oil_water_ss.count_b(solver_slau, flow)
        oil_water_ss.count_c(solver_slau, flow)
    coeff_mat = write(solver_slau.coefficient_matrix)

    s = 1
    for cell in cell_container.get_cells():
        oil_water_ss.count_a_ss(solver_slau, cell, delta_list_k[s])

    coeff_mat = write(solver_slau.coefficient_matrix)
    nevyaz_mat = write_nevyaz(solver_slau.nevyaz_vector)

    delta_list_k = sc.bicg(coeff_mat, nevyaz_mat)[0]
    #solver_slau.solve_thomas_matrix_method(0.0,0.0,0.0,0.0)
    #delta_list_k = solver_slau.get_result()

    solver_slau.clear_result()

    # for i in range(1, cell_container.get_len() - 1):
    #     cell = cell_container.get_cells()[i]
    #     cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n_plus().get_pressure_water() + delta_list_k[i - 1][0][0])
    #     cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n_plus().get_pressure_oil() + delta_list_k[i - 1][1][0])

    #update everything
    m = 1
    for i in range(0, len(delta_list_k), 2):
        cell = cell_container.get_cells()[m]
        cell.get_cell_state_n_plus().set_pressure_water(cell.get_cell_state_n_plus().get_pressure_water() + delta_list_k[i])
        m += 1

    m = 1
    for i in range(1, len(delta_list_k), 2):
        cell = cell_container.get_cells()[m]
        cell.get_cell_state_n_plus().set_pressure_oil(cell.get_cell_state_n_plus().get_pressure_oil() + delta_list_k[i])
        m += 1


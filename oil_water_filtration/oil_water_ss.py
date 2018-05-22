from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
from oil_water_filtration.OilWaterSS import OilWaterSS
from oil_water_filtration.Solver import SolverSlau

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
    delta_list_k = oil_water_ss.generate_delta_k(layer)
    for cell in cell_container.get_cells():
        cell.get_cell_state_n().set_equals_to(cell.get_cell_state_n_plus())

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
        oil_water_ss.update_pressure(cell_container, delta_list_k)

    oil_water_ss.recount_properties(cell_container)
    oil_water_ss.update_saturation(cell_container, flow_array)

    print("COUNTER " + str(counter))

    if counter == 125:
        oil_water_ss.show_results(counter, layer, cell_container)

    counter += 1
    time += oil_water_ss.tau


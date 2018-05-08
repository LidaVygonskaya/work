from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
import numpy as np
import scipy.sparse.linalg as sc
import pandas as pd
from oil_water_filtration.Enums import Components
from oil_water_filtration.OilWaterSS import OilWaterSS
from oil_water_filtration.Solver import SolverSlau
from oil_water_filtration.Write import write, write_nevyaz

np.set_printoptions(threshold=np.inf)
layer = Layer()
#solver = Solver()
oil_water_ss = OilWaterSS()
solver_slau = SolverSlau(layer.N - 2)

cell_container = CellContainer(layer.N, layer)
cell_container.initialize_cells()

flow_array = []
for i in range(cell_container.get_len() - 1):
    cells = cell_container.get_cells()
    flow_array.append(Flow(cells[i], cells[i + 1]))

for flow in flow_array:
    flow.count_cells_flows()


c = np.matrix([[flow_array[0].t_oil_water[Components.WATER.value], 0], [0, flow_array[0].t_oil_water[Components.OIL.value]]])
#print(c)
b = np.matrix([[flow_array[1].t_oil_water[Components.WATER.value], 0], [0, flow_array[1].t_oil_water[Components.OIL.value]]])
#print(b)


solver_slau.coefficient_matrix = np.zeros((98, 98, 2, 2), float)
solver_slau.nevyaz_vector = np.zeros((98, 2, 1))
for flow in flow_array:
    oil_water_ss.count_b(solver_slau, flow)
    oil_water_ss.count_c(solver_slau, flow)

for cell in cell_container.get_cells():
    oil_water_ss.count_a_ss(solver_slau, cell)


# with open('outfile.txt') as f:
#     for line in solver_slau.coefficient_matrix:
#         f.write(str(line))
coeff_mat = write(solver_slau.coefficient_matrix)
nevyaz_mat = write_nevyaz(solver_slau.nevyaz_vector)
x = sc.bicg(coeff_mat, nevyaz_mat)
#solver_slau.solve_thomas_matrix_method(0.0, 0.0, 0.0, 0.0)
solver_slau.set_zero()
print("meow")




#Water pressure


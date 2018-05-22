from oil_water_filtration.CellContainer import CellContainer
from oil_water_filtration.Flow import Flow
from oil_water_filtration.Layer import Layer
from oil_water_filtration.OilWaterIMPES import OilWaterIMPES
from oil_water_filtration.OilWaterSS import OilWaterSS
from oil_water_filtration.Solver import SolverSlau

layer = Layer()
oil_water_impes = OilWaterIMPES()
oil_water_ss = OilWaterSS()
solver_slau = SolverSlau(layer.N - 2)

cell_container = CellContainer(layer.N, layer)
cell_container.initialize_cells()

flow_array = []
for i in range(cell_container.get_len() - 1):
    cells = cell_container.get_cells()
    flow_array.append(Flow(cells[i], cells[i + 1]))

time = oil_water_impes.tau
counter = 1

import numpy as np


def write(matrix):
    toup_row = []
    for i in range(matrix.shape[0]):
        toup_col = []
        for j in range(matrix.shape[1]):
            toup_col.append(matrix[i,j])
        line_mat = np.concatenate(toup_col, axis=1)
        toup_row.append(line_mat)
    m = np.concatenate(toup_row, axis=0)
    np.savetxt('outfile_matrix.txt', m)
    return m

def write_nevyaz(matrix):
    toup_row = []
    for mat in matrix:
        toup_row.append(mat)
    m = np.concatenate(toup_row, axis=0)
    return m
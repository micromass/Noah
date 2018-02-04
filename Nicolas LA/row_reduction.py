"""Gives user the row reduced version of a matrix"""
import numpy as np

def row_reduce(matrix):
    m = len(matrix)
    n = len(matrix[0])

    for i in range(m):
        if matrix[i][i] != 0 :
            divisor = matrix[i][i]
            for j in range(n):
                matrix[i][j] = matrix[i][j]/divisor


        for row in range(m):
            cancel_num = -1 * matrix[row][i] # number which can be multiplied to cancel the term required in the matrix
            if row != i:
                matrix[row] = cancel_num * matrix[i] + matrix[row]

    print(matrix)


matrix = np.array([[1, 0, -1 ,2  , 0],
                   [-2, 1, 1, -4, 0 ],
                   [0, -1 , 2, 1, 0 ],
                   [1 , 1, 0,  4, 0 ]], float)

row_reduce(matrix)

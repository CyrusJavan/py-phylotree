"""
Utility functions
Author: Cyrus Javan
"""


def normalize(matrix):
    return_mat = matrix.copy()
    dim = len(matrix)
    total = sum(sum(v) for v in matrix)
    for i in range(dim):
        for j in range(dim):
            return_mat[i][j] = matrix[i][j] / total
    return return_mat


def normalize_list(matrix_list):
    return [normalize(mat) for mat in matrix_list]


def read_distance_matrix(file_name):
    """
    This function will parse a text-based distance matrix in the Phylip format
    param file_name: the name of the file to parse
    return: the matrix as a list of lists [[],[],[],...]
    """
    matrix_file = open(file_name, "r")
    print("Reading {}".format(file_name))
    # First line is the dimension of the NxN matrix
    dimension = int(matrix_file.readline())
    node_names = []
    return_matrix = []
    for line in matrix_file:
        tokens = line.split()
        return_matrix.append([float(s) for s in tokens[1:]])
        node_names.append(str(tokens[0]))
    return node_names, return_matrix


def combine_distance_matrices(node_name_list, matrices_list, norm=False):
    """
    This function will combine the given matrices by averaging their values
    Matrices must all be same size
    """
    num_matrices = len(matrices_list)
    dim_matrices = len(matrices_list[0])
    if norm:
        matrices_list = normalize_list(matrices_list)
    avg_matrix = [[0 for _ in range(dim_matrices)] for _ in range(dim_matrices)]
    # First sum all the correct values
    sum_matrix = [[0 for _ in range(dim_matrices)] for _ in range(dim_matrices)]
    sum_names = node_name_list[0]
    for i in range(num_matrices):
        curr_mat = matrices_list[i]
        curr_names = node_name_list[i]
        for j in range(dim_matrices):
            for k in range(dim_matrices):
                first_name = sum_names[j]
                second_name = sum_names[k]
                cj = None
                ck = None
                for l in range(dim_matrices):
                    if first_name == curr_names[l]:
                        cj = l
                    if second_name == curr_names[l]:
                        ck = l
                sum_matrix[j][k] += curr_mat[cj][ck]
    for i in range(dim_matrices):
        for j in range(dim_matrices):
            avg_matrix[i][j] = sum_matrix[i][j] / num_matrices
    return sum_names, avg_matrix


def write_distance_matrix(matrix, names, file_name):
    file = open(file_name, "w")
    file.write(str(len(matrix)) + "\n")
    for i in range(len(matrix)):
        file.write(names[i] + " ")
        for j in range(len(matrix)):
            file.write(str(matrix[i][j]) + " ")
        if i + 1 < len(matrix):
            file.write("\n")
    file.close()

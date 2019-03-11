"""
Main functionality
Author: Cyrus Javan
"""
import sys
import Utility
import ete3


def create_nodes(node_names):
    """ Returns a list of ete3.TreeNodes from the list of names """
    nodes = []
    for name in node_names:
        new_node = ete3.TreeNode(name=name)
        nodes.append(new_node)
    return nodes


def avg_dist(index, exclude_index, distance_matrix):
    """ Also called r() in the book
        Average distance of the sequence at index from all the other sequences
        except the sequence at exclude_index
    """
    running_sum = 0.0
    for i in range(len(distance_matrix)):
        if i != index and i != exclude_index:
            running_sum += distance_matrix[index][i]
    return running_sum / (len(distance_matrix) - 2)


def create_D_matrix(distance_matrix):
    """ D(i,j) = d(i,j) - (r(i) + r(j))
        r(i) = average distance to all other nodes besides i and j
    """
    dimension = len(distance_matrix)
    D_matrix = [[0 for _ in range(dimension)] for _ in range(dimension)]
    for i in range(dimension):
        for j in range(dimension):
            if i != j:
                D_matrix[i][j] = distance_matrix[i][j] - (avg_dist(i, j, distance_matrix)
                                                          + avg_dist(j, i, distance_matrix))
    return D_matrix


def add_row_col(dist_mat, dist_list):
    """
    Adds a row and column to dist_mat then returns the new matrix
    For example dist_mat =  [[0,4,5],      dist_list = [1, 2, 7]
                             [4,0,3],
                             [5,3,0]]

                Output =    [[0,4,5,1],
                             [4,0,3,2],
                             [5,3,0,7],
                             [1,2,7,0]]
    """
    return_mat = dist_mat.copy()
    for i in range(len(return_mat)):
        return_mat[i].append(dist_list[i])
    new_row = dist_list.copy()
    new_row.append(0.0)
    return_mat.append(new_row)
    return return_mat


def remove_row_col(dist_mat, indexA, indexB):
    return_mat = [[0 for _ in range(len(dist_mat) - 2)] for _ in range(len(dist_mat) - 2)]
    ii = 0
    jj = 0
    for i in range(len(dist_mat)):
        for j in range(len(dist_mat)):
            if i != indexA and i != indexB and j != indexA and j != indexB:
                return_mat[ii][jj] = dist_mat[i][j]
                jj += 1
                if jj == len(return_mat):
                    ii += 1
                    jj = 0
    return return_mat


def create_combined_dist_list(dist_mat, indexA, indexB):
    return_list = []
    for m in range(len(dist_mat)):
        if m != indexA and m != indexB:
            # new_dist formula from book
            new_dist = (dist_mat[indexA][m] + dist_mat[indexB][m] - dist_mat[indexA][indexB]) / 2.0
            return_list.append(new_dist)
    return return_list


def get_edge_length(dist_mat, indexA, indexB):
    return (dist_mat[indexA][indexB] + avg_dist(indexA, indexB, dist_mat) - avg_dist(indexB, indexA, dist_mat)) / 2.0


def create_nj_tree(node_names, distance_matrix):
    """ Create a tree using the neighbor joining algorithm
        Algorithm from: "Biological sequence analysis: Probabilistic models
        of proteins and nucleic acids" by Durbin, Eddy, Krogh and Mitchinson (1998)
        Return the root of the tree
    """
    # Create the initial list of leaf nodes
    # Our nodes will be instances of ete3.TreeNode
    # We are using ete3 so we can easily visualize the tree
    t_nodes = create_nodes(node_names)
    l_nodes = t_nodes.copy()
    dist_mat = distance_matrix.copy()

    # While there is still more than 2 leaf nodes remaining
    while len(l_nodes) > 2:
        # print("create_nj_tree::  length of leaf nodes = {}".format(len(l_nodes)))
        # "Pick a pair i, j in L for which D(i,j) is minimal"
        D_matrix = create_D_matrix(dist_mat)
        # We created the D matrix, now find the minimal value
        minimum = sys.float_info.max
        minimum_coord = (0, 0)
        for i in range(len(dist_mat)):
            for j in range(len(dist_mat)):
                if D_matrix[i][j] < minimum:
                    minimum = D_matrix[i][j]
                    minimum_coord = (i, j)
        # "Define a new node k and set d(k,m) = 1/2 * (d(i,m) + d(j,m) - d(i,j)), for all m in L except i, j"
        # We want a new distance_matrix with i, j removed and k (which is just ij) added
        # First remove i and j
        new_dist_mat = remove_row_col(dist_mat, minimum_coord[0], minimum_coord[1])
        # print("new dist mat  after remove rowCol= {}".format(new_dist_mat))
        # Add in the new node
        new_dist_list = create_combined_dist_list(dist_mat, minimum_coord[0], minimum_coord[1])
        # print("new dist list = {}".format(new_dist_list))
        new_dist_mat = add_row_col(new_dist_mat, new_dist_list)
        # print("new dist mat  after add RowCol = {}".format(new_dist_mat))
        # So now we have a new_dist_mat with k added, and i and j removed
        # Now, create a TreeNode with i and j as children and set the correct distances
        k = ete3.TreeNode()
        dist_i_k = get_edge_length(dist_mat, minimum_coord[0], minimum_coord[1])
        dist_j_k = dist_mat[minimum_coord[0]][minimum_coord[1]] - dist_i_k
        k.add_child(child=l_nodes[minimum_coord[0]], dist=dist_i_k)
        k.add_child(child=l_nodes[minimum_coord[1]], dist=dist_j_k)

        # Remove i and j, add k
        del l_nodes[minimum_coord[0]]
        if minimum_coord[0] > minimum_coord[1]:
            del l_nodes[minimum_coord[1]]
        else:
            del l_nodes[minimum_coord[1] - 1]
        l_nodes.append(k)

        # print("old dist mat = {}".format(dist_mat))
        # print("new dist mat = {}".format(new_dist_mat))
        # print("length old dist mat = {}".format(len(dist_mat[0])))
        # print("length new dist mat = {}".format(len(new_dist_mat[0])))
        # finally update the dist_mat
        dist_mat = new_dist_mat
    # Now there is only 2 leaves remaining
    # Define a new root
    root = ete3.TreeNode(dist=0)
    root.add_child(child=l_nodes[0], dist=dist_mat[0][1]/2)
    root.add_child(child=l_nodes[1], dist=dist_mat[0][1] / 2)
    # Choose either of the nodes to be the root
    # root = l_nodes[0]
    # root.add_child(child=l_nodes[1], dist=dist_mat[0][1])
    return root


def layout(node):
    if node.is_leaf():
        N = ete3.AttrFace("name", fsize=30)
        ete3.faces.add_face_to_node(N, node, 0, position="aligned")


def _main():
    print("Reading matrix from file")
    root = None
    output_file_name = None
    if len(sys.argv) < 2:
        print("Please supply at least one path to a distance matrix")
        quit(1)
    elif sys.argv[1] == "matrix":  # Usage: python3 PhyloTree.py matrix OUTPUT_FILE_NAME [PATH_TO_INPUT_MATRICES . . .]
        node_names_list = []
        dist_mat_list = []
        for mat in sys.argv[3:]:
            (node_names, dist_mat) = Utility.read_distance_matrix(mat)
            node_names_list.append(node_names)
            dist_mat_list.append(dist_mat)
        (node_names, dist_mat) = Utility.combine_distance_matrices(node_names_list, dist_mat_list, norm=False)
        print("Writing combined matrix to {}".format(sys.argv[2]))
        Utility.write_distance_matrix(dist_mat, node_names, sys.argv[2])
        exit(0)
    elif len(sys.argv) == 2:  # Only one matrix provided
        # Read the distance matrix in to memory
        (node_names, dist_mat) = Utility.read_distance_matrix(sys.argv[1])
        print("Creating tree with Neighbor-Joining from 1 matrix")
        # Our neighbor joining algorithm creates rooted trees
        # the root is created a node in the middle of the final edge added to the tree
        root = create_nj_tree(node_names, dist_mat)
        gene_name = sys.argv[1].split("/")[-1]
        output_file_name = "{}_tree.png".format(gene_name)
    else:  # more than one matrix, average them then make tree
        node_names_list = []
        dist_mat_list = []
        for mat in sys.argv[1:]:
            (node_names, dist_mat) = Utility.read_distance_matrix(mat)
            node_names_list.append(node_names)
            dist_mat_list.append(dist_mat)
        (node_names, dist_mat) = Utility.combine_distance_matrices(node_names_list, dist_mat_list, norm=False)
        print("Creating tree with Neighbor-Joining from {} matrices".format(len(sys.argv[1:])))
        root = create_nj_tree(node_names, dist_mat)
        system_name = sys.argv[1].split("/")[-1]
        output_file_name = "{}_combined_tree.png".format(system_name)

    ts = ete3.TreeStyle()
    ts.show_branch_length = True
    ts.title.add_face(ete3.TextFace(output_file_name[0:-4], fsize=20), column=0)
    nstyle = ete3.NodeStyle()
    nstyle["shape"] = "sphere"
    nstyle["size"] = 6
    nstyle["fgcolor"] = "darkred"
    for n in root.traverse():
        n.set_style(nstyle)
    root.show(tree_style=ts)
    root.render(output_file_name, tree_style=ts, units="px", w=600)


if __name__ == '__main__':
    _main()

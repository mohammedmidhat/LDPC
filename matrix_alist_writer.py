## Writes a matrix in the alist format into a text file
##
##Author: Mohammed Al Ai Baky
##Created: 9/24/2017



import numpy as np


def matrix_to_alist(H, file_name):
    f = open(file_name, 'w')

    r = int(H.shape[0])
    c = int(H.shape[1])

    f.write(str(c)+" "+str(r)+"\n")

    col_weights_arr = np.count_nonzero(H, axis=0)
    row_weights_arr = np.count_nonzero(H, axis=1)

    col_weight_max = max(col_weights_arr)
    row_weight_max = max(row_weights_arr)

    f.write(str(col_weight_max)+" "+str(row_weight_max)+"\n");
    for i in col_weights_arr:
        f.write(str(i)+ " ")
    f.write('\n')
    for i in row_weights_arr:
        f.write(str(i)+ " ")
    f.write('\n')
    
    for i in range(c):
        indices = np.where(H[:,i]==1)[0]
        for j in indices:
            f.write(str(j+1)+" ")
        for j in range(col_weight_max-col_weights_arr[i]):
            f.write("0 ")
        f.write("\n")

    for i in range(r):
        indices = np.where(H[i,:]==1)[0]
        for j in indices:
            f.write(str(j+1)+" ")
        for j in range(row_weight_max-row_weights_arr[i]):
            f.write("0 ")
        f.write("\n")

    f.close()

def alist_to_mat(filename):
    f = open(filename)

    line = f.readline()
    line_content = line.split(' ')
    c = int(line_content[0])
    r = int(line_content[1])
    H = np.array([[0 for j in range(c)] for i in range(r)])

    line = f.readline()
    line_content = line.split(' ')
    col_weight_max = int(line_content[0])
    row_weight_max = int(line_content[1])

    line = f.readline()
    col_weight_arr = line.split(' ')
    col_weight_arr = col_weight_arr[:-1]

    line = f.readline()
    row_weight_arr = line.split(' ')
    row_weight_arr = row_weight_arr[:-1]

    for i in range(c):
        line = f.readline()
        indices = line.split(' ')
        indices = indices[:-1]
        for j in range(int(col_weight_arr[i])):
            H[int(indices[j])-1,i] = 1
    
    f.close()

    return H

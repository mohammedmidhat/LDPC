import array
import numpy as np

H_file = open('204.33.484 (N=204,K=102,M=102,R=0.5).txt')
H_string = H_file.read()


def make_numpy_H_columns(string, k, n):
    array_result = array.array('l',[[0 for j in range(n)] for i in range(n-k)])
    #array_result = [[0 for j in range(n)] for i in range(n-k)]
    list_columns = string.split('\n')
    for i in range(len(list_columns)):
        for j in list_columns[i].split('\t'):
            array_result[int(j)-1][i] = 1
    return array_result


def make_numpy_H_rows(string, k, n):
    array_result = array.array('l',[[0 for j in range(n)] for i in range(n-k)])
    #array_result = [[0 for j in range(n)] for i in range(n-k)]
    list_columns = string.split('\n')
    for i in range(len(list_columns)):
        print('i = ' + str(i))
        for j in list_columns[i].split('\t'):
            print('j = ' + str(j))
            array_result[i][int(j)-1] = 1
    return array_result


def check_weights_type(string, k, n):
    list_weights = string.split('\n')
    if(len(list_weights) == k):
        print 'rows'
    elif(len(list_weights) == n):
        print 'columns'
    return

def get_102_204_3_6():
    return make_numpy_H_columns(H_string, 102, 204)

def exp():
    return [[1,2],[3,4]]

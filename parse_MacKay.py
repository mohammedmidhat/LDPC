import numpy as np

####H_file = open('204.33.484 (N=204,K=102,M=102,R=0.5).txt')
H_file = open('2.txt')
H_string = H_file.read()

##another_file = open('1.txt')
##another_string = another_file.read()

##def make_numpy_H_columns(string, k, n):
##    array = np.array([[0 for j in range(n)] for i in range(n-k)])
####    array = [[0 for j in range(n)] for i in range(n-k)]
##    list_columns = string.split('\n')
##    for i in range(len(list_columns)):
##        for j in list_columns[i].split('\t'):
##            array[int(j)-1][i] = 1
##    return array
##
##
def make_numpy_H_rows(string, k, n):
    array = np.array([[0 for j in range(n)] for i in range(n-k)])
##    array = [[0 for j in range(n)] for i in range(n-k)]
    list_columns = string.split('\n')
    for i in range(len(list_columns)):
        print('i = ' + str(i))
        for j in list_columns[i].split('\t'):
            print('j = ' + str(j))
            array[i][int(j)-1] = 1
    return array


def check_weights_type(string, k, n):
    list_weights = string.split('\n')
    if(len(list_weights) == k):
        print 'rows'
    elif(len(list_weights) == n):
        print 'columns'
    return


def exp():
    #return 'yep'
    return np.array([1,2])


import numpy as np
import copy


def find_smallest(array):
    if len(array) == 1:
        return 0
    elif len(array) == 2:
        if array[0] <= array[1]:
            return 0
        else:
            return 1
    else:
        arrayA = array[:len(array)/2]
        arrayB = array[(len(array)/2):]
        smallA = find_smallest(arrayA)
        smallB = find_smallest(arrayB)
        if arrayA[smallA] <= arrayB[smallB]:
            return smallA
        else:
            return len(arrayA) + smallB


class peg():

    """
    Progressive edge growth algorithm for generating
    LDPC matrices. The algorithm is obtained from [1]
    """

    def __init__(self, nvar, nchk, degree_sequence):
        self.degree_sequence = degree_sequence
        self.nvar = nvar
        self.nchk = nchk
        self.H = np.zeros((nchk, nvar), dtype = np.int32)
        self.chk_degrees = np.zeros(nchk, dtype = np.int32)
        self.girth = 1000
    
    def grow_edge(self, var, chk):
        self.H[chk, var] = 1
        self.chk_degrees[chk] += 1

    def bfs(self, var):
        var_list = np.zeros(self.nvar, dtype = np.int32)
        var_list[var] = 1
        cur_chk_list = []
        new_chk_list = []
        for i in range(self.nchk):
            cur_chk_list.append(0)
            new_chk_list.append(0)
        chk_Q = []
        var_Q = []
        var_Q.append(var)
        l = 0
        while(True):
            for _vars in var_Q:
                for i in range(self.nchk):
                    if self.H[i, _vars] == 1:
                        if cur_chk_list[i] == 0:
                            new_chk_list[i] = 1
                            chk_Q.append(i)
            var_Q = []
            for _chks in chk_Q:
                for j in range(self.nvar):
                    if self.H[_chks, j] == 1:
                        if var_list[j] == 0:
                            var_list[j] = 1
                            var_Q.append(j)
            chk_Q = []
            if new_chk_list.count(1) == self.nchk:
                new_chk = self.find_smallest_chk(cur_chk_list)
                self.girth = l+2 if l+2 < self.girth else self.girth
                return new_chk
            elif np.array_equal(new_chk_list, cur_chk_list):
                new_chk = self.find_smallest_chk(cur_chk_list)
                return new_chk
            else:
                cur_chk_list = copy.copy(new_chk_list)

    def find_smallest_chk(self, cur_chk_list):
        index = []
        degree = []
        for i in range(len(cur_chk_list)):
            if cur_chk_list[i] == 0:
                index.append(i)
                degree.append(self.chk_degrees[i])
        return index[find_smallest(degree)]

    def progressive_edge_growth(self):
        for var in range(self.nvar):
            print "edge growth at var", var
            for k in range(self.degree_sequence[var]):
                if k == 0:
                    smallest_degree_chk = find_smallest(self.chk_degrees)
                    self.grow_edge(var, smallest_degree_chk)
                else:
                    chk = self.bfs(var)
                    self.grow_edge(var, chk)


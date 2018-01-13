import numpy as np
import copy


class peg():

    """
    Progressive edge growth algorithm for generating
    LDPC matrices. The algorithm is obtained from [1]
    """

    def __init__(self, nvar, nchk, var_degree_sequence, chk_degree_sequence, p):
        self.var_degree_sequence = var_degree_sequence
        self.chk_degree_sequence = chk_degree_sequence
        self.cur_var_degree_sequence = np.zeros(p*nvar, dtype = np.int32)
        self.cur_chk_degree_sequence = np.zeros(p*nchk, dtype = np.int32)
        self.nvar = nvar
        self.nchk = nchk
        self.p = p
        self.H = np.zeros((p*nchk, p*nvar), dtype = np.int32)
        self.girth = 1000

    # returns the node of the smallest degree
    # returns one at random if there are more than one
    def find_nodes_of_smallest_degree(self, array):
        sorted_array = np.argsort(array)
        for i in sorted_array:
            if(array[i] < self.chk_degree_sequence[i/self.p]):
                return i

    def bfs(self, var):
        var_list = np.zeros(self.nvar*self.p, dtype = np.int32)
        var_list[var] = 1
        cur_chk_list = np.array([0 for i in range(self.nchk*self.p)])
        new_chk_list = np.array([0 for i in range(self.nchk*self.p)])
        chk_Q = []
        var_Q = []
        var_Q.append(var)
        l = 0
        while(True):
            for _vars in var_Q:
                neighbors = np.where(self.H[:,_vars]==1)[0]
                for i in neighbors:
                    if cur_chk_list[i] == 0:
                        new_chk_list[i] = 1
                        chk_Q.append(i)
            var_Q = []
            for _chks in chk_Q:
                neighbors = np.where(self.H[_chks,:]==1)[0]
                for j in neighbors:
                    if var_list[j] == 0:
                        var_list[j] = 1
                        var_Q.append(j)
            chk_Q = []
            
            if np.count_nonzero(new_chk_list) == self.nchk*self.p:
                new_chk = self.find_smallest_chk(cur_chk_list)
                self.girth = l+2 if l+2 < self.girth else self.girth
                return new_chk
            elif np.array_equal(new_chk_list, cur_chk_list):
                new_chk = self.find_smallest_chk(cur_chk_list)
                return new_chk
            else:
                cur_chk_list = copy.copy(new_chk_list)

    def find_smallest_chk(self, array):
        ## Mark all the entire group of each node in the tree
        input_array_len = len(array)
        input_array = copy.copy(array)
        for i in range(0,input_array_len,self.p):
            for j in range(i, i+self.p):
                if(input_array[j]!=0):
                    for k in range(i, i+self.p):
                        input_array[k] = -1
                    break
        
        
        valid_chks = np.where(input_array==0)[0]
        chk_degree_sequence_expanded = np.array([j for j in self.chk_degree_sequence for i in range(self.p)])
        sorted_array = np.argsort(chk_degree_sequence_expanded)
        for i in sorted_array:
            if((i*self.p in valid_chks) and (array[i] < self.chk_degree_sequence[i/self.p])):
                return i

    def grow_first_edge_from_group(self, var):
        _chk = self.find_nodes_of_smallest_degree(self.cur_chk_degree_sequence)
        self.H[_chk,var] = 1
        self.cur_chk_degree_sequence[_chk] += 1
        self.cur_var_degree_sequence[var] += 1

        for i in range(1,self.p):
            chk_node = (_chk/self.p)*self.p + (_chk+i)%self.p
            self.H[chk_node, var+i] = 1
            self.cur_chk_degree_sequence[chk_node] += 1
            self.cur_var_degree_sequence[var+i] += 1

    def grow_edge_from_group(self, var, chk):
        self.H[chk,var] = 1
        self.cur_chk_degree_sequence[chk] += 1
        self.cur_var_degree_sequence[var] += 1

        for i in range(1,self.p):
            chk_node = (chk/self.p)*self.p + (chk+i)%self.p
            self.H[chk_node, var+i] = 1
            self.cur_chk_degree_sequence[chk_node] += 1
            self.cur_var_degree_sequence[var+i] += 1

    def progressive_edge_growth(self):
        for var in range(self.nvar):
            print "edge growth at var", var
            for k in range(self.var_degree_sequence[var]):
                if k == 0:
                    self.grow_first_edge_from_group(var*self.p)
                else:
                    edge_growth_result = self.bfs(var*self.p)
                    print edge_growth_result
                    self.grow_edge_from_group(var*self.p, edge_growth_result)

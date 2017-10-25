import numpy as np
from scipy import sparse


protograph_1 = {0:[4,2,3],1:[0,1,2],2:[1,2,3],3:[1,3,4],4:[0,1,4],5:[0,4,5],6:[4,3,5],7:[0,4,1],8:[0,5,3],9:[0,2,5],10:[2,5,3],11:[5,2,1]}


def make_H_b(protograph,n,k,num_checks):
    H = np.array([[0 for i in range(n)] for i in range(n-k)])
    num_copies = n/len(protograph)
    #print('num_copies = ' + str(num_copies))
    protograph_copies = {}
    for i in protograph:
        num_edges = len(protograph[i])
        #print('num_edges = ' + str(num_edges))
        perm = [np.random.permutation(num_copies) for k in range(num_edges)]
        #print('perm_'+str(i)+' = '+str(perm))
        for j in range(num_copies):
            protograph_copies[j*len(protograph)+i] = [num_checks*perm[k][j]+protograph[i][k] for k in range(num_edges)]
            #print protograph_copies
    #print protograph_copies
    for i in protograph_copies:
        H[protograph_copies[i],i] = 1
    return H


def reduce_redundant_H(H):
    redundant_rows_indices = []
    for i in range(len(H)-1):
        if(i not in redundant_rows_indices):
            for j in range(i+1,len(H)):
                if(j not in redundant_rows_indices):
                    if(np.all(H[i]==H[j])):
                        redundant_rows_indices.append(j)
    return redundant_rows_indices


def make_H_nb(protograph,n,k,q):
    H_nb = ''
    return H_nb

mat = make_H_b(protograph_1,204,102,6).astype('int')
np.savetxt('proto_based.txt',mat,delimiter=',',newline='\n')
